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

ellip_corr_path = os.path.join(dirname(dirname(dirname(abspath(__file__)))), 'ellip-corr')
sys.path.append(ellip_corr_path)
from PyEllipCorr import PyEllipCorr
from seismic.pick_harvester.travel_time import TTInterpolator
import traceback

class SSSTRelocator(ParametricData):
    def __init__(self, csv_catalog, auto_pick_files=[], auto_pick_phases=[],
                 events_only=False, phase_list='P Pg Pb Pn S Sg Sb Sn', temp_dir='./'):
        super(SSSTRelocator, self).__init__(csv_catalog, auto_pick_files, auto_pick_phases,
                                            events_only, phase_list, temp_dir)
        self.EARTH_RADIUS = 6371. #km
        self.DEG2KM = np.pi/180 * self.EARTH_RADIUS
        self.kdtree = None
        self._lonlatalt2xyz = None
        self.ecdist = None
        self.edepth_km = None
        self.ecolat = None
        self.azim = None
        self.elev_corr = None
        self.ellip_corr = None
        self.vsurf = None
        self._tti = TTInterpolator()
        self._pyellipcorr = PyEllipCorr()
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

    def compute_essentials(self):
        elon = self.events['lon'][self.event_id_to_idx[self.arrivals['event_id']]]
        elat = self.events['lat'][self.event_id_to_idx[self.arrivals['event_id']]]
        slon = self.arrivals['lon']
        slat = self.arrivals['lat']

        azim, _, ecdist = self._geod.inv(elon, elat, slon, slat)
        self.azim = np.array(azim.astype('f4'))
        self.ecdist = np.array(ecdist.astype('f4'))

        self.edepth_km = self.events['depth_km'][self.event_id_to_idx[self.arrivals['event_id']]]
        self.ecolat = 90 - self.events['lat'][self.event_id_to_idx[self.arrivals['event_id']]]
    # end func

    def compute_elevation_correction(self):
        vsurf = self.vsurf
        elev_km = self.arrivals['elev_m'] / 1e3

        dtdd = self._tti.get_dtdd(self.arrivals['phase'],
                                  self.ecdist,
                                  self.edepth_km)
        elev_corr = vsurf * (dtdd / self.DEG2KM)
        elev_corr = np.power(elev_corr, 2.)
        elev_corr[elev_corr > 1.] = 1./elev_corr[elev_corr > 1.]
        elev_corr = np.sqrt(1. - elev_corr)
        elev_corr = elev_corr * elev_km / vsurf

        elev_corr = np.array(elev_corr.astype('f4'))
        elev_corr[np.isnan(elev_corr)] = 0.
        self.elev_corr = elev_corr
    # end func

    def compute_ellipticity_correction(self):
        ellip_corr = np.zeros(len(self.arrivals), dtype='f4')

        for i, arrival in enumerate(self.arrivals):
            ellip_corr[i] = self._pyellipcorr.get_correction(arrival['phase'], self.ecdist[i],
                                                             self.edepth_km[i], self.ecolat[i],
                                                             self.azim[i])
        # end for
        ellip_corr[np.isnan(ellip_corr)] = 0.
        self.ellip_corr = ellip_corr
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
        global_values = np.empty(nelem, dtype=dtype)

        fn = self._temp_dir + '/{}.npy'.format(self.rank)
        np.save(fn, rank_values)
        self.comm.Barrier()

        for irank in np.arange(self.nproc):
            b, e = displacements[irank], displacements[irank]+counts[irank]
            fn = self._temp_dir + '/{}.npy'.format(irank)
            global_values[b:e] = np.load(fn.format(irank))
        # end for
        self.comm.Barrier()

        return global_values
    # end func
# end class

if __name__ == "__main__":
    if(0):
        sr = SSSTRelocator('./small_merge_catalogues_output.csv',
                           auto_pick_files=['small_p_combined.txt', 'small_s_combined.txt'],
                           auto_pick_phases=['P', 'S'],
                           events_only=False)
        print('computing distance, azimuth, elev_corr')
        sr.compute_essentials()
        sr.compute_elevation_correction()
        sr.compute_ellipticity_correction()
        if(sr.rank == 0): print(sr.azim, sr.ecdist, sr.elev_corr, sr.ellip_corr)
        test_ellip_corr = sr._sync_var(split_list(sr.ellip_corr, sr.nproc)[sr.rank])
        assert np.all(sr.ellip_corr == test_ellip_corr), 'Sync-var failed'
    else:
        sr = SSSTRelocator('./merge_catalogues_output.csv',
                           auto_pick_files=['p_combined.txt', 's_combined.txt'],
                           auto_pick_phases=['P', 'S'],
                           events_only=False)
        print('computing distance, azimuth, elev_corr')
        sr.compute_essentials()
        sr.compute_elevation_correction()
        sr.compute_ellipticity_correction()
        if(sr.rank == 0): print(sr.azim, sr.ecdist, sr.elev_corr, sr.ellip_corr)
        test_ellip_corr = sr._sync_var(split_list(sr.ellip_corr, sr.nproc)[sr.rank])
        assert np.all(sr.ellip_corr == test_ellip_corr), 'Sync-var failed'
    # end if
# end if
