from mpi4py import MPI
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
from os.path import abspath, dirname

ellip_corr_path = os.path.join(dirname(dirname(dirname(abspath(__file__)))), 'ellip-corr')
sys.path.append(ellip_corr_path)
from PyEllipCorr import PyEllipCorr
from seismic.pick_harvester.travel_time import TTInterpolator

class SSSTRelocator(ParametricData):
    def __init__(self, csv_catalog, auto_pick_files=[], auto_pick_phases=[],
                 events_only=False, phase_list='P Pg Pb P* Pn S Sg Sb S* Sn'):
        super(SSSTRelocator, self).__init__(csv_catalog, auto_pick_files, auto_pick_phases,
                                            events_only, phase_list)

        self.kdtree = None
        self._lonlatalt2xyz = None
        self.local_events_indices = np.array(split_list(np.arange(len(self.events)), self.nproc)[self.rank], dtype=np.int)
        self.local_arrivals_indices = np.array(split_list(np.arange(len(self.arrivals)), self.nproc)[self.rank], dtype=np.int)

        self._sync_events()
        self._sync_arrivals()
        self._initialize_kdtree()

        self.tti = TTInterpolator()
        self.ellipcorr = PyEllipCorr()
    # end func

    def _initialize_kdtree(self, ellipsoidal_distance=False):
        ER = 6371e3 #m

        elons = self.events['lon']
        elats = self.events['lat']
        ealts = -self.events['depth_km'] * 1e3
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

    def _sync_events(self):
        # sync events across ranks
        event_counts = np.array(self.comm.allgather(len(self.local_events_indices)))

        nelem = np.sum(event_counts)
        displacements = np.zeros(self.nproc)
        displacements[1:] = np.cumsum(event_counts[:-1])
        global_events = np.empty(nelem, dtype=self.event_fields)
        type_map = {'i4': MPI.INT, 'f4': MPI.FLOAT, 'f8': MPI.DOUBLE}
        for name, dtype in zip(self.event_fields['names'], self.event_fields['formats']):
            if(dtype in ['i4', 'f4', 'f8']):
                temp = np.zeros(nelem, dtype=dtype)
                self.comm.Allgatherv(np.array(self.events[name][self.local_events_indices]),
                                     [temp, event_counts,
                                      displacements, type_map[dtype]])
                global_events[name][:] = temp
            else:
                length = np.dtype(dtype).itemsize
                temp = np.empty(nelem*length, dtype='b')
                self.comm.Allgatherv(np.array(self.events[name][self.local_events_indices]).tobytes(),
                                     [temp, event_counts*length,
                                      displacements*length, MPI.CHAR])
                global_events[name] = np.frombuffer(temp, dtype=dtype)
            # end if
        # end for
        self.events = global_events
    # end func

    def _sync_arrivals(self):
        # sync arrivals across all ranks
        arrival_counts = np.array(self.comm.allgather(len(self.local_arrivals_indices)))

        nelem = np.sum(arrival_counts)
        displacements = np.zeros(self.nproc)
        displacements[1:] = np.cumsum(arrival_counts[:-1])
        global_arrivals = np.empty(nelem, dtype=self.arrival_fields)
        type_map = {'i4': MPI.INT, 'f4': MPI.FLOAT, 'f8': MPI.DOUBLE}
        for name, dtype in zip(self.arrival_fields['names'], self.arrival_fields['formats']):
            if(dtype in ['i4', 'f4', 'f8']):
                temp = np.zeros(nelem, dtype=dtype)
                self.comm.Allgatherv(np.array(self.arrivals[name][self.local_arrivals_indices]),
                                     [temp, arrival_counts,
                                      displacements, type_map[dtype]])
                global_arrivals[name][:] = temp
            else:
                length = np.dtype(dtype).itemsize
                temp = np.empty(nelem*length, dtype='b')
                self.comm.Allgatherv(np.array(self.arrivals[name][self.local_arrivals_indices]).tobytes(),
                                     [temp, arrival_counts*length,
                                      displacements*length, MPI.CHAR])
                global_arrivals[name] = np.frombuffer(temp, dtype=dtype)
            # end if
        # end for
        self.arrivals = global_arrivals
    # end func
# end class

if __name__ == "__main__":
    if(1):
        sr = SSSTRelocator('./small_merge_catalogues_output.csv',
                           auto_pick_files=['small_p_combined.txt', 'small_s_combined.txt'],
                           auto_pick_phases=['P', 'S'],
                           events_only=False)
    else:
        sr = SSSTRelocator('./merge_catalogues_output.csv',
                           auto_pick_files=['old_p_combined.txt', 'old_s_combined.txt'],
                           auto_pick_phases=['P', 'S'],
                           events_only=False)
    # end if
# end if