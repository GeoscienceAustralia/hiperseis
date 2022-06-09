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
from mpi4py import MPI

class SSSTRelocator(ParametricData):
    def __init__(self, csv_catalog, auto_pick_files=[], auto_pick_phases=[],
                 events_only=False, phase_list='P Pg Pb P* Pn S Sg Sb S* Sn'):
        super(SSSTRelocator, self).__init__(csv_catalog, auto_pick_files, auto_pick_phases,
                                            events_only, phase_list)
        self.EARTH_RADIUS = 6371 #km
        self.DEG2KM = np.pi/180 * self.EARTH_RADIUS
        self.kdtree = None
        self._lonlatalt2xyz = None
        self.local_events_indices = np.array(split_list(np.arange(len(self.events)), self.nproc)[self.rank], dtype='i4')
        self.local_arrivals_indices = np.array(split_list(np.arange(len(self.arrivals)), self.nproc)[self.rank], dtype='i4')
        self.ecdist = None
        self.azim = None
        self.elev_corr = None
        self.vsurf = None
        self._tti = TTInterpolator()
        self._ellipcorr = PyEllipCorr()
        self._geod = Geod(a=180/np.pi, f=0)

        self._initialize_kdtree()

        # initialize surface velocity array
        vs_surf = 3.46        # Sg velocity km/s for elevation corrections
        vp_surf = 5.8         # Pg velocity km/s for elevation corrections
        pidx = self.is_p()
        self.vsurf = np.zeros(len(self.arrivals), dtype='f4')
        self.vsurf[pidx] = vp_surf
        self.vsurf[~pidx] = vs_surf
    # end func

    def is_p(self):
        result = self.arrivals['phase'].astype('S1') == b'P'
        return result
    # end func

    def compute_distance_azimuth(self):
        elon = self.events['lon'][self.event_id_to_idx[self.arrivals['event_id'][self.local_arrivals_indices]]]
        elat = self.events['lat'][self.event_id_to_idx[self.arrivals['event_id'][self.local_arrivals_indices]]]
        slon = self.arrivals['lon'][self.local_arrivals_indices]
        slat = self.arrivals['lat'][self.local_arrivals_indices]

        local_azim, _, local_ecdist = self._geod.inv(elon, elat, slon, slat)
        self.azim = self._sync_var(np.array(local_azim.astype('f4')))
        self.ecdist = self._sync_var(np.array(local_ecdist.astype('f4')))
    # end func

    def compute_elevation_correction(self):
        vsurf = self.vsurf[self.local_arrivals_indices]
        print(self.rank, 'hereA')
        #elev_km = (self.arrivals['elev_m'][self.local_arrivals_indices]) / 1e3
        elev_km = np.ones(len(vsurf), dtype='f4')
        print(self.rank, 'hereB')
        print(self.rank, 'here0')
        #edepth_km = self.events['depth_km'][self.event_id_to_idx[self.arrivals['event_id'][self.local_arrivals_indices]]]
        edepth_km = np.ones(len(vsurf), dtype='f4')*100
        print(self.rank, 'here1')
        dtdd = self._tti.get_dtdd(self.arrivals['phase'][self.local_arrivals_indices],
                                  self.ecdist[self.local_arrivals_indices],
                                  edepth_km)

        elev_corr = vsurf * (dtdd / self.DEG2KM)
        elev_corr *= elev_corr
        elev_corr[elev_corr > 1.] = 1./elev_corr[elev_corr > 1.]
        elev_corr = np.sqrt(1. - elev_corr)
        elev_corr *= elev_km / vsurf
        print(self.rank, 'here2')

        local_elev_corr = np.array(elev_corr.astype('f4'))
        print(self.rank, 'here3', np.sum(np.isnan(local_elev_corr)))
        self.elev_corr = self._sync_var(local_elev_corr)
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
# end class

if __name__ == "__main__":
    if(0):
        sr = SSSTRelocator('./small_merge_catalogues_output.csv',
                           auto_pick_files=['small_p_combined.txt', 'small_s_combined.txt'],
                           auto_pick_phases=['P', 'S'],
                           events_only=False)
        if(sr.rank==0): print('computing distance, azimuth, elev_corr')
        sr.compute_distance_azimuth()
        sr.compute_elevation_correction()
        if(sr.rank==0): print(sr.azim, sr.ecdist, sr.elev_corr)
    else:
        sr = SSSTRelocator('./merge_catalogues_output.csv',
                           auto_pick_files=['p_combined.txt', 's_combined.txt'],
                           auto_pick_phases=['P', 'S'],
                           events_only=False)
        if(sr.rank==0):
            print(len(sr.arrivals), sr.arrivals.dtype.itemsize)
            print('computing distance, azimuth, elev_corr')
        sr.compute_distance_azimuth()
        sr.compute_elevation_correction()
        if(sr.rank==0): print(sr.azim, sr.ecdist, sr.elev_corr)
    # end if
# end if