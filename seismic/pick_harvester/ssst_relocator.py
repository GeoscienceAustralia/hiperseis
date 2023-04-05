"""
Description:
    Implements the SSSTRelocator class that relocates events and redefines arrival phases

References:

CreationDate:   20/07/22
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     20/07/22   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from ordered_set import OrderedSet as set
import numpy as np
from scipy.spatial import cKDTree
import pyproj
from collections import defaultdict
from seismic.pick_harvester.utils import split_list

from seismic.pick_harvester.parametric_data import ParametricData, EARTH_RADIUS_KM
import os
from tqdm import tqdm
from seismic.pick_harvester.travel_time import TTInterpolator
from seismic.pick_harvester.ellipticity import Ellipticity
from seismic.pick_harvester.ssst_utils import compute_ellipticity_correction, \
                                                compute_elevation_correction, \
                                                VP_SURF, VS_SURF
import h5py
import json

class SSSTRelocator(ParametricData):
    def __init__(self, csv_catalog, auto_pick_files=[], auto_pick_phases=[], config_fn=None,
                 phase_list='P Pg Pb Pn S Sg Sb Sn',
                 p_residual_cutoff=5, s_residual_cutoff=10, temp_dir='./'):
        """
        :param csv_catalog: path to catalog file in csv format
        :param auto_pick_files: optional list of files containing automatic arrivals. Each file
                                should contain arrivals of the same phase (e.g. P and S arrivals
                                output from pick.py)
        :param auto_pick_phases: optional list of phases corresponding to each file provided in
                                 auto_pick_files
        :param config_fn: a config file as described in ../Readme.md
        :param phase_list: a space-separated list of phases to be read from the catalogue -- all
                           other phases not in this list are discarded
        :param p_residual_cutoff: ± residual cutoff for P-arrivals -- P-phases with residuals
                                  exceeding this cutoff are discarded.
        :param s_residual_cutoff: ± residual cutoff for S-arrivals -- S-phases with residuals
                                  exceeding this cutoff are discarded
        :param temp_dir: path to a temporary folder to be used for syncing data across processors.
                         Note that this temporary folder must be accessible by all MPI ranks, e.g.
                         within a project folder on the NCI.
        """
        super(SSSTRelocator, self).__init__(csv_catalog, auto_pick_files, auto_pick_phases,
                                            False, phase_list, temp_dir)
        self.config_fn = config_fn
        self.p_residual_cutoff = p_residual_cutoff
        self.s_residual_cutoff = s_residual_cutoff
        self._lonlatalt2xyz = None
        self.vsurf = None
        self.residual = np.zeros(len(self.arrivals), dtype='f4')
        self._source_enum = defaultdict(int)
        self._arrival_source = None
        self._tti = TTInterpolator()
        self._ellipticity = Ellipticity()
        self.station_id = None
        self.tcorr = np.zeros(len(self.arrivals), dtype='f4')
        self.sst_tcorr = np.zeros(len(self.arrivals), dtype='f4')
        self.event_quality = np.ones(len(self.events), dtype='?') # all events are considered good at the outset

        self.events_imask = np.zeros(len(self.events), dtype='?') # events to be relocated
        self.arrivals_imask = np.zeros(len(self.arrivals), dtype='?') # arrivals whose phases are to be redefined
        self.arrivals_relocated = np.zeros(len(self.arrivals), dtype='?') # arrivals that belong to events being relocated

        # keep a copy of original events for comparing relocations against
        self.orig_events = np.array(self.events)

        # keep a copy of input arrival phases for producing relevant plots
        self.orig_phases = np.array(self.arrivals['phase'])

        # a copy of residuals computed before phases are redefined and ssst-relocations applied
        self.r0 = None

        self._assign_station_ids()
        self._label_data()
        self._initialize_spatial_functors()

        # parse config file and lable events/arrivals to be processed
        self._parse_config()

        # initialize surface velocity array (only needs to be computed once, since
        # elevations are fixed and P/S arrivals are not interchanged)
        pidx = self.is_P(self.arrivals['phase'])
        sidx = self.is_S(self.arrivals['phase'])
        self.vsurf = np.zeros(len(self.arrivals), dtype='f4')
        self.vsurf[pidx] = VP_SURF
        self.vsurf[sidx] = VS_SURF
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

        @param phase: dim(self.local_arrival_indices)
        @param imask: dim(self.local_arrival_indices)
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
        tcorr = self.tcorr[self.local_arrivals_indices][imask]

        azim, _, ecdist = self._geod.inv(elon, elat, slon, slat)

        elev_corr = compute_elevation_correction(self._tti, phase, vsurf, elev_km, ecdist, edepth_km)
        ellip_corr = compute_ellipticity_correction(self._ellipticity, phase, ecdist, edepth_km, elat, azim)
        ptt = self._tti.get_tt(phase, ecdist, edepth_km)

        # compute empirical tt
        ett = self.arrivals['arrival_ts'][self.local_arrivals_indices][imask] - \
              self.events['origin_ts'][indices][imask]

        residual = ett - tcorr - (ptt.astype('f8') + elev_corr + ellip_corr)

        # set invalid residuals
        residual[ptt == self._tti.fill_value] = self._tti.fill_value

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

        @param imask: dim(self.local_arrival_indices)
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

        local_new_phase[p_indices[np.fabs(local_new_residual[p_indices]) > self.p_residual_cutoff]] = b'Px'
        local_new_phase[s_indices[np.fabs(local_new_residual[s_indices]) > self.s_residual_cutoff]] = b'Sx'

        local_new_residual[p_indices[np.fabs(local_new_residual[p_indices]) > self.p_residual_cutoff]] = self._tti.fill_value
        local_new_residual[s_indices[np.fabs(local_new_residual[s_indices]) > self.s_residual_cutoff]] = self._tti.fill_value

        self.arrivals['phase'] = self._sync_var(local_new_phase)
        self.residual = self._sync_var(local_new_residual)
    # end func

    def compute_sst_correction(self, min_slope_ratio=5):
        """
        Operates on global data
        @param min_slope_ratio:
        """
        imask = self.arrivals_imask

        stas = np.unique(self.station_id)
        local_sids = split_list(stas, self.nproc)[self.rank]
        local_corrs = np.ones(len(self.arrivals), dtype='f4') * np.nan

        for sid in tqdm(local_sids, desc='SST Rank {}: '.format(self.rank)):
            sta_arrival_indices = np.argwhere((self.station_id == sid) & (imask)).flatten()

            # drop automatic arrival indices where slope_ratio < min_slope_ratio
            sta_arrival_indices = sta_arrival_indices[ ((self.is_AUTO_arrival[sta_arrival_indices]) & \
                                                        (self.arrivals['quality_measure_slope'][sta_arrival_indices] >= min_slope_ratio)) | \
                                                       (~self.is_AUTO_arrival[sta_arrival_indices]) ]

            sta_arrivals = self.arrivals[sta_arrival_indices]

            phases = set(sta_arrivals['phase'])
            for p in phases:
                pimask = sta_arrivals['phase'] == p

                residuals = self.residual[sta_arrival_indices][pimask]
                residuals = np.ma.masked_array(residuals, mask=residuals==self._tti.fill_value)
                sst_corr = 0
                if(np.sum(~residuals.mask)): sst_corr = np.ma.mean(residuals)

                local_corrs[sta_arrival_indices[pimask]] = sst_corr
            # end for
        # end for

        indices = np.argwhere(~np.isnan(local_corrs)).flatten()
        self.sst_tcorr = self._gather(len(self.arrivals), local_corrs[indices], indices)
        self.tcorr[:] = self.sst_tcorr[:]
    # end func

    def compute_ssst_correction(self, ball_radius_km=55, min_slope_ratio=5):
        """
        Operates on global data
        @param ball_radius_km: radius of sphere for computing SSST corrections
        """
        imask = self.arrivals_imask

        stas = np.unique(self.station_id)
        local_sids = split_list(stas, self.nproc)[self.rank]
        local_corrs = np.ones(len(self.arrivals), dtype='f4') * np.nan

        MIN_SSST_ARRIVALS = 5 # ssst corrections are computed when number of arrivals >= this value
        for sid in tqdm(local_sids, desc='SSST Rank {}: '.format(self.rank)):
            sta_arrival_indices = np.argwhere((self.station_id == sid) & (imask)).flatten()

            # drop automatic arrival indices where slope_ratio < min_slope_ratio
            sta_arrival_indices = sta_arrival_indices[ ((self.is_AUTO_arrival[sta_arrival_indices]) & \
                                                        (self.arrivals['quality_measure_slope'][sta_arrival_indices] >= min_slope_ratio)) | \
                                                       (~self.is_AUTO_arrival[sta_arrival_indices]) ]

            sta_arrivals = self.arrivals[sta_arrival_indices]

            phases = set(sta_arrivals['phase'])
            phase_trees = {}
            phase_residuals = {}
            phase_arrival_indices = {}
            for p in phases:
                pimask = sta_arrivals['phase'] == p

                event_indices = self.event_id_to_idx[sta_arrivals['event_id'][pimask]]

                # drop events where event_quality is False
                event_indices = event_indices[self.event_quality[event_indices]]

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

                if(np.sum(~nbr_phase_residuals.mask) >= MIN_SSST_ARRIVALS):
                    corr = np.ma.median(nbr_phase_residuals)
                    local_corrs[sta_arrival_indices[i]] = corr
                # end if
            # end for
        # end for

        indices = np.argwhere(~np.isnan(local_corrs)).flatten()
        global_tcorr = self._gather(len(self.arrivals), local_corrs[indices], indices)

        # ===========================================================
        # Time corrections (both sst and ssst) applied to arrivals
        # that belong to events being relocated are cumulative, since
        # they result in relocation of respective event origins.
        # Arrivals whose phases are simply redefined without their
        # respective events being relocated, on the other hand, must
        # factor in both sst and ssst time corrections
        # for the current ssst-relocation iteration.
        # ===========================================================
        self.tcorr[self.arrivals_relocated] = global_tcorr[self.arrivals_relocated]
        self.tcorr[~self.arrivals_relocated] = self.sst_tcorr[~self.arrivals_relocated] + \
                                               global_tcorr[~self.arrivals_relocated]
    # end func

    def relocate_events(self, imask=None, reloc_niter=10, reloc_dlon=0.25,
                        reloc_dlat=0.25, reloc_ddep=10, reloc_sfrac=0.5,
                        min_slope_ratio=5):

        """
        Operates on local data

        @param imask: dim(self.local_event_indices)
        @param reloc_niter:
        @param reloc_dlon:
        @param reloc_dlat:
        @param reloc_ddep:
        @param reloc_sfrac:
        @return:
        """
        SOL_DLON = 0.25 # +/- 0.25 degrees
        SOL_DLAT = 0.25 # +/- 0.25 degrees
        SOL_DDEPTH_KM = 20 # +/- 20 km
        SOL_DT = 5 # +/- 5 seconds
        MIN_ARRIVALS = 3 # minimum number of usable arrivals an event must have to be considered of good quality

        if (imask is None): imask = np.ones(len(self.local_events_indices), dtype='?')

        qual = np.ones(len(self.local_events_indices), dtype='?')
        selected_local_indices = np.argwhere(imask).flatten()
        for i, iev in enumerate(tqdm(self.local_events_indices[imask], desc='RELOC Rank {}: '.format(self.rank))):
            event = self.events[iev]

            arrival_indices = np.argwhere(self.arrivals['event_id'] == event['event_id']).flatten()
            #b=len(arrival_indices)
            # drop automatic arrival indices where slope_ratio < min_slope_ratio
            arrival_indices = arrival_indices[ ((self.is_AUTO_arrival[arrival_indices]) & \
                                                (self.arrivals['quality_measure_slope'][arrival_indices] >= min_slope_ratio)) | \
                                               (~self.is_AUTO_arrival[arrival_indices]) ]
            #a=len(arrival_indices)
            #print('{}/{}'.format(a,b))

            ones = np.ones(len(arrival_indices))

            slon = self.arrivals['lon'][arrival_indices]
            slat = self.arrivals['lat'][arrival_indices]
            vsurf = self.vsurf[arrival_indices]
            elev_km = self.arrivals['elev_m'][arrival_indices] / 1e3
            phase = self.arrivals['phase'][arrival_indices]
            ett = self.arrivals['arrival_ts'][arrival_indices] - event['origin_ts']
            tcorr = self.tcorr[arrival_indices]

            elon0 = event['lon']
            elat0 = event['lat']
            edepth_km0 = event['depth_km']

            elon_best = elon0
            elat_best = elat0
            edepth_km_best = edepth_km0

            dlon = reloc_dlon
            dlat = reloc_dlat
            ddep = reloc_ddep
            sfrac = reloc_sfrac

            dt = 0
            dtbest = dt
            npick_used_best = 0

            E0 = 0
            Ebest = 1e20
            for ireloc in np.arange(reloc_niter):
                for ilon in np.arange(-1, 2):
                    elon = elon0 + ilon * dlon
                    for ilat in np.arange(-1, 2):
                        elat = elat0 + ilat * dlat
                        for idep in np.arange(-1, 2):
                            edepth_km = edepth_km0 + idep * ddep
                            azim, _, ecdist = self._geod.inv(elon * ones, elat * ones, slon, slat)

                            elev_corr = compute_elevation_correction(self._tti, phase, vsurf,
                                                                     elev_km, ecdist, edepth_km * ones)
                            ellip_corr = compute_ellipticity_correction(self._ellipticity, phase, ecdist,
                                                                        edepth_km * ones, elat * ones, azim)
                            tt = self._tti.get_tt(phase, ecdist, edepth_km * ones)
                            tt = np.ma.masked_array(tt, mask=tt == self._tti.fill_value)

                            npick_used = np.sum(~tt.mask)
                            if (npick_used):
                                ptt = tt + elev_corr + ellip_corr
                                residual = ett - tcorr - ptt

                                residmed = np.ma.mean(residual)
                                E0 = np.sqrt(np.mean(np.power(residual - residmed, 2)))

                                if (E0 < Ebest):
                                    Ebest = E0
                                    elon_best = elon
                                    elat_best = elat
                                    edepth_km_best = edepth_km
                                    dtbest = residmed
                                    npick_used_best = npick_used
                                # end if
                            # end if
                        # end for
                    # end for
                # end for
                elon0, elat0, edepth_km0 = elon_best, elat_best, edepth_km_best
                dlon, dlat, ddep = np.array([dlon, dlat, ddep]) * sfrac
            # end for
            self.events[iev]['lon'] = elon_best
            self.events[iev]['lat'] = elat_best
            self.events[iev]['depth_km'] = edepth_km_best
            self.events[iev]['origin_ts'] += dtbest

            good = (np.fabs(self.orig_events[iev]['lon'] - self.events[iev]['lon']) < SOL_DLON) & \
                   (np.fabs(self.orig_events[iev]['lat'] - self.events[iev]['lat']) < SOL_DLAT) & \
                   (np.fabs(self.orig_events[iev]['depth_km'] - self.events[iev]['depth_km']) < SOL_DDEPTH_KM) & \
                   ((np.fabs(self.orig_events[iev]['origin_ts'] - self.events[iev]['origin_ts']) < SOL_DT))
            qual[selected_local_indices[i]] = (npick_used_best >= MIN_ARRIVALS) & (good)
        # end for

        self.events['lon'] = self._sync_var(self.events['lon'][self.local_events_indices])
        self.events['lat'] = self._sync_var(self.events['lat'][self.local_events_indices])
        self.events['depth_km'] = self._sync_var(self.events['depth_km'][self.local_events_indices])
        self.events['origin_ts'] = self._sync_var(self.events['origin_ts'][self.local_events_indices])

        self.event_quality = self._sync_var(qual)
    # end func

    def ssst_relocate(self, ssst_niter=5, ball_radius_km=55, min_slope_ratio=5, output_fn=None):
        """
        Operates both on global and local data
        @param ssst_niter:
        @param ball_radius_km:
        @param output_fn:
        @return:
        """

        def dump_h5(fn, iter):
            h = h5py.File(fn, 'a')

            eg = h.create_group('{}/events'.format(iter))
            ag = h.create_group('{}/arrivals'.format(iter))

            # dump events
            for var in self.event_fields['names']:
                eg.create_dataset(var, data=self.events[var])
            # end for

            # dump event-quality
            eg.create_dataset('event_quality', data=self.event_quality)

            # dump arrivals
            for var in self.arrival_fields['names']:
                ag.create_dataset(var, data=self.arrivals[var])
            # end for

            # dump ssst-tcorrs
            else: ag.create_dataset('tcorr', data=self.tcorr)

            # dump residual
            ag.create_dataset('residual', data=self.residual)

            if(iter == 0):
                # dump r0
                ag.create_dataset('r0', data=self.r0)

                # dump all is_SRC_[event/arrival] boolean arrays -- note that these flags
                # remain immutable since creation
                vlist = [v for v in vars(self) if 'is' in v and ('event' in v or 'arrival' in v)]
                for v in vlist:
                    if('event' in v):
                        eg.create_dataset(v, data=getattr(self, v))
                    elif('arrival' in v):
                        ag.create_dataset(v, data=getattr(self, v))
                    # end if
                # end for

                # dump event-id-to-idx
                eg.create_dataset('event_id_to_idx', data=self.event_id_to_idx)
            # end if
            h.close()
        # end func

        # =======================================================
        # compute all residuals and save r0, since some
        # residuals are recomputed in redefine_phases
        # =======================================================
        self.compute_residual()
        self.r0 = np.array(self.residual)

        if(1):
            #===========================================================
            # SST relocate events and redefine phases
            #===========================================================
            self.redefine_phases(imask=self.arrivals_imask[self.local_arrivals_indices])
            self.compute_sst_correction(min_slope_ratio=min_slope_ratio)
            self.relocate_events(imask=self.events_imask[self.local_events_indices],
                                 min_slope_ratio=min_slope_ratio)

            # dump results
            if (self.rank == 0 and output_fn): dump_h5(output_fn, 0)
        # end if

        if(1):
            # ===========================================================
            # SSST relocate events and redefine phases
            # ===========================================================
            for issst in tqdm(np.arange(1, ssst_niter+1), desc='SSST-iter: '):
                self.compute_residual()
                self.redefine_phases(imask=self.arrivals_imask[self.local_arrivals_indices])

                self.compute_ssst_correction(ball_radius_km=ball_radius_km,
                                             min_slope_ratio=min_slope_ratio)

                self.relocate_events(imask=self.events_imask[self.local_events_indices],
                                     min_slope_ratio=min_slope_ratio)

                if(issst == ssst_niter):
                    self.compute_residual()
                    self.redefine_phases(imask=self.arrivals_imask[self.local_arrivals_indices])
                # end if

                # dump results
                if (self.rank == 0 and output_fn): dump_h5(output_fn, issst)
            # end for
        # end if
    # end func

    def _initialize_spatial_functors(self, ellipsoidal_distance=False):
        if(self.rank == 0): print('Initializing spatial functors..')

        ER = EARTH_RADIUS_KM * 1e3 #m

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

    def _gather(self, final_size, rank_values, rank_indices):
        assert(len(rank_values) == len(rank_indices))

        dtype = rank_values.dtype
        global_values = np.zeros(final_size, dtype=dtype)

        fn = self._temp_dir + '/gather.h5'

        for irank in np.arange(self.nproc):
            if(self.rank == irank):
                hf = h5py.File(fn, 'a')
                dset = hf.create_dataset("{}".format(self.rank), rank_values.shape, dtype=rank_values.dtype)
                dset[:] = rank_values
                dset_i = hf.create_dataset("{}_i".format(self.rank), rank_values.shape, dtype=rank_indices.dtype)
                dset_i[:] = rank_indices
                hf.close()
            # end if
            self.comm.Barrier()
        # end for

        hf = h5py.File(fn, 'r')
        for irank in np.arange(self.nproc):
            dset = hf['{}'.format(irank)][:]
            dset_i = hf['{}_i'.format(irank)][:]
            global_values[dset_i] = dset
        # end for
        hf.close()

        self.comm.Barrier()
        if(self.rank == 0): os.remove(fn)

        return global_values
    # end func

    def _parse_config(self):
        if(self.rank == 0): print('Parsing config file: {}..'.format(self.config_fn))

        c = json.load(open(self.config_fn))

        sources = set([i.decode() for i in set(self.events['source'])])

        # sanity check 1: ensure keys are available
        if (len(set(c.keys()) - sources)):
            assert 0, 'Keys {} not found in data set. Aborting..'.format(list(set(c.keys()) - sources))
        # end if

        # sanity check 2: ensure config blocks are specified for each source
        if (len(sources - set(c.keys()))):
            assert 0, 'Config block not found for source: {}. Aborting..'.format(list(sources - set(c.keys())))
        # end if

        # sanity check 3: each block must have three entries (events, preexisting_arrivals, automatic_arrivals)
        items = set(['events', 'preexisting_arrivals', 'automatic_arrivals'])
        for k in c.keys():
            if (set(c[k].keys()) != items):
                assert 0, 'Invalid keys found: {}'.format(list(items.symmetric_difference(set(c[k].keys()))))
            # end if
        # end for

        # sanity check 4: check entries for events, preexisting_arrivals and automatic_arrivals
        expected_values = {'events': ['relocate', 'fixed'],
                           'preexisting_arrivals': ['redefine', 'fixed'],
                           'automatic_arrivals': ['redefine', 'fixed']}
        for k in c.keys():
            for ek, ev in expected_values.items():
                if (c[k][ek] not in ev):
                    assert 0, 'Invalid values found for "{}" in block "{}". Expected values are: {}'.format(ek, k, ev)
                # end if
            # end for
        # end for

        for k in c.keys():
            event_attr = 'is_{}_event'.format(k)
            arrival_attr = 'is_{}_arrival'.format(k)

            if (c[k]['events'] == 'relocate'):
                self.events_imask |= getattr(self, event_attr)
            # end if

            if (c[k]['preexisting_arrivals'] == 'redefine'):
                self.arrivals_imask |= ((getattr(self, arrival_attr)) & (~self.is_AUTO_arrival))
            # end if

            if (c[k]['automatic_arrivals'] == 'redefine'):
                self.arrivals_imask |= ((getattr(self, arrival_attr)) & (self.is_AUTO_arrival))
            # end if
        # end for

        # imask for arrivals that belong to events being relocated
        self.arrivals_relocated = np.isin(self.arrivals['event_id'], self.events['event_id'][self.events_imask])
    # end func
# end class
