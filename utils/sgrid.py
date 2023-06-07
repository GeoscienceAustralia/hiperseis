"""
Description:
   Implements a generic SGrid reader (binary/ascii)

CreationDate:   15/02/22
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     15/02/22   RH

"""
import os
import numpy as np
from scipy.spatial import cKDTree
from pyproj import Proj
from seismic.misc import rtp2xyz

class SGrid:
    def __init__(self, filename, utm_zone):
        self._fn = filename
        self._utm_zone = utm_zone
        self._is_binary = False
        self._prop_file = None
        self._points_file = None
        self._prop_alignment = None
        self._nx = None
        self._ny = None
        self._nz = None
        self._zscale = 1.
        self._data = None
        self._tree = None
        self._proj = None
        self._EARTH_RADIUS = 6371e3

        self._setup_projection()
        self._read_header()
        self._read_data()
    # end func

    def get_values(self, lon, lat, depth, nnbr=1, p=4):
        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)
        d = np.atleast_1d(depth)

        rtp = np.zeros((lon.shape[0], 3))
        rtp[:,0] = self._EARTH_RADIUS - d
        rtp[:,1] = np.radians(lat - 90)
        rtp[:,2] = np.radians(lon)

        xyz = rtp2xyz(rtp[:,0], rtp[:,1], rtp[:,2])

        ds, ids = self._tree.query(xyz, k=nnbr)

        if(nnbr>1):
            idw = 1./np.power(ds, p)
            vals = np.sum(idw * np.reshape(self._data[ids.flatten(), 3], ids.shape), axis=1) / \
                   np.sum(idw, axis=1)
        else:
            vals = self._data[ids, 3]
        # end if

        return vals
    # end func

    def __str__(self):
        val = """file-name: {}\n
                 is-binary: {}\n
                 points-file: {}\n
                 prop-file: {}\n
                 nx: {}\n
                 ny: {}\n
                 nz: {}\n
                 zscale: {}
                 projection: {}
                 """.format(self._fn, self._is_binary, self._points_file,
                            self._prop_file, self._nx, self._ny, self._nz,
                            self._zscale, self._proj.to_wkt())
        return val
    # end func

    def _read_header(self):
        for line in open(self._fn).readlines():
            if(line.startswith('ascii')):
                if(line.split(':')[1].strip().lower() == 'on'):
                    self._is_binary = False
                else: self._is_binary = True
            # end if

            if(self._is_binary):
                if(line.startswith('PROP_FILE')): self._prop_file = line.split(' ')[2].strip()
                if(line.startswith('POINTS_FILE')): self._points_file = line.split(' ')[1].strip()
            else:
                if(line.startswith('ASCII_DATA_FILE')): self._prop_file = line.split(' ')[1].strip()
            # end if

            if(line.startswith('AXIS_N ')):
                self._nx, self._ny, self._nz = map(int, line.strip().split(' ')[1:])
            # end if

            if(line.startswith('PROP_ALIGNMENT')):
                self._prop_alignment = line.split(' ')[-1].strip()
            # end if

            if(line.startswith('ZPOSITIVE')):
                # convert elevation to depth
                if(line.split(' ')[-1].strip().lower() == 'elevation'): self._zscale *= -1
            # end if
        # end for

        self._prop_file = os.path.join(os.path.dirname(self._fn), self._prop_file)
        if(self._points_file): self._points_file = os.path.join(os.path.dirname(self._fn),
                                                                self._points_file)
    # end func

    def _read_data(self):
        if(self._is_binary):
            points = np.fromfile(self._points_file, dtype='>f4')
            vals = np.fromfile(self._prop_file, dtype='>f4')

            points = points.reshape(self._nz, self._ny, self._nx, 3)
            if(self._prop_alignment == 'CELLS'):
                cpoints = np.zeros((self._nz-1, self._ny-1, self._nx-1, 3))
                cpoints[:,:,:,0] = (points[:-1, :-1, :-1, 0] + points[1:, 1:, 1:, 0])/2.
                cpoints[:,:,:,1] = (points[:-1, :-1, :-1, 1] + points[1:, 1:, 1:, 1])/2.
                cpoints[:,:,:,2] = (points[:-1, :-1, :-1, 2] + points[1:, 1:, 1:, 2])/2.

                points = cpoints
            # end if

            points = points.reshape((-1, 3))

            assert points.shape[0] == vals.shape[0], \
            'Array size mismatch detected between coordinates and values in {}'. format(self._fn)

            self._data = np.hstack((points, vals[:,None]))
        else:
            self._data = np.loadtxt(self._prop_file, comments='*')
            assert self._data.shape[0] == self._nx * self._ny * self._nz, \
            'Error loading data from: {}'.format(self._prop_file)
        # end if

        self._data = self._data[:,:4]

        self._data[:,0], self._data[:,1] = self._proj(self._data[:,0], self._data[:,1], inverse=True)
        self._data[:,2] *= self._zscale

        rtp = np.zeros((self._data.shape[0], 3))
        rtp[:,0] = self._EARTH_RADIUS - self._data[:,2]
        rtp[:,1] = np.radians(self._data[:,1] - 90)
        rtp[:,2] = np.radians(self._data[:,0])
        xyz = rtp2xyz(rtp[:, 0], rtp[:,1], rtp[:,2])
        self._tree = cKDTree(xyz, balanced_tree=False)
    # end func

    def _setup_projection(self):
        # set up projection
        try:
            north = None
            if(self._utm_zone[-1].lower() == 'n'): north = True
            elif(self._utm_zone[-1].lower() == 's'): north = False
            else: raise ValueError('Zone number must be followed by n/s, example 53s')

            zone = int(self._utm_zone[:-1])
            self._proj = Proj('+proj=utm +zone={} +{} +datum=WGS84 +units=m +no_defs'.format(
                               zone, 'north' if north else 'south'))
        except Exception as e:
            print(e)
            print('UTM zone should be provided as INTn/s, e.g. 53s')
        # end try
    # end func
# end class
