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
        self._lonlatalt2xyz = None
        self.vsurf = None
        self.residual = None
        self._source_enum = defaultdict(int)
        self._arrival_source = None
        self._tti = TTInterpolator()
        self._ellipticity = Ellipticity()
        self.station_id = None
        self.ssst_tcorr = np.zeros(len(self.arrivals), dtype='f4')

        self._assign_station_ids()
        self._label_data()
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
        netsta_ids = self.comm.bcast(netsta_ids, root=0)

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
        ssst_tcorr = self.ssst_tcorr[self.local_arrivals_indices][imask]

        azim, _, ecdist = self._geod.inv(elon, elat, slon, slat)

        elev_corr = self.compute_elevation_correction(phase, vsurf, elev_km, ecdist, edepth_km)
        ellip_corr = self.compute_ellipticity_correction(phase, ecdist, edepth_km, elat, azim)
        ptt = self._tti.get_tt(phase, ecdist, edepth_km)

        # compute empirical tt
        ett = self.arrivals['arrival_ts'][self.local_arrivals_indices][imask] - \
              self.events['origin_ts'][indices][imask]

        residual = ett - ssst_tcorr - (ptt.astype('f8') + elev_corr + ellip_corr)
        #residual = ett - (ptt.astype('f8') + elev_corr + ellip_corr)

        # set invalid residuals
        residual[ptt == self._tti.fill_value] = self._tti.fill_value
        #print("Rank: {}, Num invalid residuals: {}".format(self.rank, np.sum(ptt==self._tti.fill_value)))

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

        local_new_phase[p_indices] = self.p_phases[np.argmin(np.fabs(p_residuals), axis=1)]
        local_new_phase[s_indices] = self.s_phases[np.argmin(np.fabs(s_residuals), axis=1)]

        p_argmin_indices = np.argmin(np.fabs(p_residuals), axis=1)
        s_argmin_indices = np.argmin(np.fabs(s_residuals), axis=1)
        local_new_residual[p_indices] = np.take_along_axis(p_residuals, p_argmin_indices[:, None], axis=1).flatten()
        local_new_residual[s_indices] = np.take_along_axis(s_residuals, s_argmin_indices[:, None], axis=1).flatten()

        p_threshold = 5 #TODO
        s_threshold = 10

        #print('total {}, bad {}'.format(len(p_indices), np.sum(np.fabs(local_new_residual[p_indices]) > p_threshold)))
        #print('total {}, bad {}'.format(len(s_indices), np.sum(np.fabs(local_new_residual[s_indices]) > s_threshold)))

        local_new_phase[p_indices[np.fabs(local_new_residual[p_indices]) > p_threshold]] = b'Px'
        local_new_phase[s_indices[np.fabs(local_new_residual[s_indices]) > s_threshold]] = b'Sx'

        local_new_residual[p_indices[np.fabs(local_new_residual[p_indices]) > p_threshold]] = self._tti.fill_value
        local_new_residual[s_indices[np.fabs(local_new_residual[s_indices]) > s_threshold]] = self._tti.fill_value

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

    def compute_ssst_correction(self, imask=None, ball_radius_km=11.):
        """
        Operates on global data
        @param imask: dim(len(arrivals))
        """
        if(imask is None): imask = np.ones(len(self.arrivals), dtype='?')

        #test_imask = (~self.is_AUTO_arrival) | ((self.is_AUTO_arrival) & (self.arrivals['quality_measure_slope']>=2))
        #imask = imask & test_imask

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
                                                                             ball_radius_km * 1e3)
                nbr_phase_residuals = phase_residuals[sta_arrival['phase']][indices]
                nbr_phase_residuals = np.ma.masked_array(nbr_phase_residuals,
                                                         mask=nbr_phase_residuals==self._tti.fill_value)

                if(np.sum(~nbr_phase_residuals.mask) >= 5): #TODO
                    corr = np.ma.median(nbr_phase_residuals)
                    local_corrs[sta_arrival_indices[i]] = corr
                # end if
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

    def _label_data(self):
        """
        Label events and arrivals
        """
        sources = set(self.events['source'])
        for i, source in enumerate(sources): self._source_enum[source] = i + 1

        if(self.rank == 0): print('Event sources found {}..'.format(self._source_enum.items()))
        if(self.rank == 0): print('Labelling events by event source..')

        # create source-type attributes for arrivals
        for source in sources:
            setattr(self, 'is_{}_event'.format(source.decode()), self.events['source'] == source)
        # end for

        if(self.rank == 0): print('Labelling arrivals by event source..')

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

        # create source-type attributes for arrivals
        for source in sources:
            setattr(self, 'is_{}_arrival'.format(source.decode()), self._arrival_source == self._source_enum[source])
        # end for

        # attribute marking automatic picks
        setattr(self, 'is_AUTO_arrival', self.arrivals['quality_measure_slope'] > -1)
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

def relocate_events(sr, imask=None):
    debug=False

    if (imask is None): imask = np.ones(len(sr.local_events_indices), dtype='?')

    reloc_niter = 10
    for iev in tqdm(sr.local_events_indices[imask]):
        event = sr.events[iev]

        #if(event['event_id'] != 29711): continue
        #if(event['event_id'] != 27161): continue

        arrival_indices = np.argwhere(sr.arrivals['event_id'] == event['event_id']).flatten()

        if debug: print(event)
        if debug: print(sr.arrivals[arrival_indices])

        if(0):
            temp = []
            for ai in arrival_indices:
                if(sr.is_AUTO_arrival[ai]):
                    pass
                    if(sr.arrivals['quality_measure_slope'][ai]>=2): temp.append(ai)
                else:
                    temp.append(ai)
                # end if
            # end for
            arrival_indices = np.array(temp, dtype='i4')
        # end if

        ones = np.ones(len(arrival_indices))

        slon = sr.arrivals['lon'][arrival_indices]
        slat = sr.arrivals['lat'][arrival_indices]
        vsurf = sr.vsurf[arrival_indices]
        elev_km = sr.arrivals['elev_m'][arrival_indices] / 1e3
        phase = sr.arrivals['phase'][arrival_indices]
        ett = sr.arrivals['arrival_ts'][arrival_indices] - event['origin_ts']
        ssst_tcorr = sr.ssst_tcorr[arrival_indices]

        elon0 = event['lon']
        elat0 = event['lat']
        edepth_km0 = event['depth_km']

        elon_best = elon0
        elat_best = elat0
        edepth_km_best = edepth_km0

        dt = 0
        dtbest = dt
        npick_used_best = len(arrival_indices)

        dlat = 0.25
        dlon = 0.25
        ddep = 10  # km
        sfrac = 0.5

        E0 = 0
        Ebest = 1e20
        ibest=None
        #print(elon_best, elat_best, edepth_km_best, dtbest, npick_used_best)
        for ireloc in np.arange(reloc_niter):
            for ilon in np.arange(-1, 2):
                elon = elon0 + ilon * dlon
                for ilat in np.arange(-1, 2):
                    elat = elat0 + ilat * dlat
                    for idep in np.arange(-1, 2):
                        edepth_km = edepth_km0 + idep * ddep
                        azim, _, ecdist = sr._geod.inv(elon * ones, elat * ones, slon, slat)

                        elev_corr = sr.compute_elevation_correction(phase, vsurf, elev_km, ecdist, edepth_km * ones)
                        ellip_corr = sr.compute_ellipticity_correction(phase, ecdist, edepth_km * ones, elat * ones,
                                                                       azim)
                        tt = sr._tti.get_tt(phase, ecdist, edepth_km * ones)
                        tt = np.ma.masked_array(tt, mask=tt == sr._tti.fill_value)

                        npick_used = np.sum(~tt.mask)
                        if (npick_used):
                            ptt = tt + elev_corr + ellip_corr
                            residual = ett - ssst_tcorr - ptt

                            residmed = np.ma.mean(residual)
                            E0 = np.sqrt(np.mean(np.power(residual - residmed, 2)))

                            if (E0 < Ebest):
                                if debug: print(E0, Ebest, npick_used,'/',len(tt), '----------')
                                Ebest = E0
                                elon_best = elon
                                elat_best = elat
                                edepth_km_best = edepth_km
                                dtbest = residmed
                                npick_used_best = npick_used
                                ibest = (ireloc, ilon, ilat, idep)
                            # end if
                        # end if
                    # end for
                # end for
            # end for
            elon0, elat0, edepth_km0 = elon_best, elat_best, edepth_km_best
            dlon, dlat, ddep = np.array([dlon, dlat, ddep]) * sfrac
        # end for
        #print(elon_best, elat_best, edepth_km_best, dtbest, npick_used_best)
        #print('-------')
        #results.append([elon_best, elat_best, edepth_km_best, dtbest, npick_used_best])
        sr.events[iev]['lon'] = elon_best
        sr.events[iev]['lat'] = elat_best
        sr.events[iev]['depth_km'] = edepth_km_best
        sr.events[iev]['origin_ts'] += dtbest
    # end for
    if debug: print(sr.events[event['event_id'] == 29711])

    sr.events['lon'] = sr._sync_var(sr.events['lon'][sr.local_events_indices])
    sr.events['lat'] = sr._sync_var(sr.events['lat'][sr.local_events_indices])
    sr.events['depth_km'] = sr._sync_var(sr.events['depth_km'][sr.local_events_indices])
    sr.events['origin_ts'] = sr._sync_var(sr.events['origin_ts'][sr.local_events_indices])
# end func

def ssst_relocate():
    ssst_niter = 5

    sr = SSSTRelocator('./merge_catalogues_output.csv',
                       auto_pick_files=['p_combined.txt', 's_combined.txt'],
                       auto_pick_phases=['P', 'S'],
                       events_only=False)

    ball_radius_km = 55
    for issst in tqdm(np.arange(ssst_niter), desc='SSST-iter: '):
        sr.compute_residual()
        sr.redefine_phases(imask=sr.is_GA_arrival[sr.local_arrivals_indices])

        sr.compute_ssst_correction(imask=sr.is_GA_arrival,
                                   ball_radius_km=float(ball_radius_km)/float(issst+1))

        np.save('ssst_tcorr_{}.npy'.format(issst), sr.ssst_tcorr)
        np.save('relocated_ga_events_{}.npy'.format(issst), sr.events)
        np.save('relocated_ga_arrivals_{}.npy'.format(issst), sr.arrivals)
        np.save('relocated_ga_residuals_{}.npy'.format(issst), sr.residual)

        relocate_events(sr, sr.is_GA_event[sr.local_events_indices])
    # end for
# end func

if __name__ == "__main__":
    ssst_relocate()
    exit(0)


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
        sr.redefine_phases(imask=sr.is_GA_arrival[sr.local_arrivals_indices] & sr.is_AUTO_arrival[sr.local_arrivals_indices])

        if(sr.rank==0): print('Computing SSST corrections..')
        sr.compute_ssst_correction()
        if(sr.rank==0):
            np.save('test_arrivals.npy', sr.arrivals)
            np.save('ssst_tcorr.npy', sr.ssst_tcorr)
        # end if

        if(sr.rank == 0):
            filt = sr.is_GA_arrival & sr.is_AUTO_arrival
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
