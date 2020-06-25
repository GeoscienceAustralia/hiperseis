# -*- coding: utf-8 -*-
"""
Class for parsing Alexei's tomographic parametrisation file as an object.

Created on Wed Apr 03 13:52:28 2019

@author: Marcus W. Haynes
"""

import numpy as np
from scipy.interpolate import griddata


class grid():
    def __init__(self, in_file='param'):
        self.file_name = in_file
        self._parse_parametrisation()

    def _parse_line(self, l):
        """
        Returns a list of white-space delimited entries, removing any
        comments indicated by an `!'.
        """

        l = l.strip('\n').strip('\r')
        if '!' in l:
            l = l.split('!')[0]
        if l == '':
            return False
        while '  ' in l:
            l = l.replace('  ',' ')
        if l[0] == ' ':
            l = l[1:]
        if l[-1] == ' ':
            l = l[:-1]
        return l.split(' ')

    def _parse_parametrisation(self):
        """
        Reads in the local and global grid parameters from an external
        parametrisation file (by default we load  param_file='./param')
        """

        with open(self.file_name, 'r') as F:
            # Rewind to file beginning
            F.seek(0)

            ## PART 1 - GLOBAL GRID
            # Read-in the global grid (number of cells and resolution)
            gnx, gny, gnz, gdx, gdy = self._parse_line(F.readline())
            self.gnx = int(gnx)
            self.gny = int(gny)
            self.gnz = int(gnz)
            self.gdx = float(gdx)
            self.gdy = float(gdy)

            # There is currently no option to reduce the coverage of the global
            # grid, so resolution information is a little redundant. Nevertheless,
            # let's perform a sanity test.
            assert (360. / self.gdx) == self.gnx
            assert (180. / self.gdy) == self.gny

            # Load global depth nodes
            self.gdepth = np.zeros(self.gnz + 1)
            for i in range(self.gnz + 1):
                try:
                    depth = np.array(self._parse_line(F.readline()), dtype=float)
                except ValueError as err:
                    print(('There appear to be less global depth levels than ' + \
                           'was indicated! The header suggested %d levels, but ' + \
                           'I can only parse %d. Please check the param file.') % (self.gnz + 1, i))
                    raise err
                self.gdepth[i] = depth

            # Global depth in meters
            self.gmeters = 1000. * self.gdepth

            ## PART 2 - LOCAL GRID
            # Read-in the local grid (extent and number of cells)
            try:
                LON_min, LON_max, LAT_min, LAT_max, nx, ny, nz = self._parse_line(F.readline())
                LON_min = float(LON_min)
                LON_max = float(LON_max)
                LAT_min = float(LAT_min)
                LAT_max = float(LAT_max)
                self.LON = np.array([LON_min, LON_max])
                self.LAT = np.array([LAT_min, LAT_max])
                self.nx = int(nx)
                self.ny = int(ny)
                self.nz = int(nz)
            except ValueError as err:
                print(('There appear to be more global depth levels than ' + \
                       'was indicated! The header suggested %d levels; ' + \
                       'please check the param file.') % (self.nz + 1))
                raise err

            # Local grid resolution
            self.dx = (self.LON[1] - self.LON[0]) / self.nx
            self.dy = (self.LAT[1] - self.LAT[0]) / self.ny

            # Load local depth nodes
            self.rdepth = np.zeros(self.nz + 1)
            for i in range(self.nz + 1):
                try:
                    depth = np.array(self._parse_line(F.readline()), dtype=float)
                except ValueError as err:
                    print(('There appear to be less local depth levels than ' + \
                           'was indicated! The header suggested %d levels, but ' + \
                           'I can only parse %d. Please check the param file.') % (self.nz + 1, i))
                    raise err
                self.rdepth[i] = depth

            # Local depth in meters
            self.rmeters = 1000. * self.rdepth

            # Get block IDs
            self.nmax, self.no_params = np.array(self._parse_line(F.readline()), dtype=int)

            # Get grid
            self.grid = np.array(self._parse_line(F.readline()), dtype=int)

    def save(self):
        """
        Write the parametrisation file to disk
        """

        with open(self.file_name, 'w') as F:
            F.write(' %11d %11d %11d %11.3f %11.3f\n' % (self.gnx, self.gny, self.gnz, self.gdx, self.gdy))
            for i in range(self.gnz+1):
                F.write(' %11.3f\n' % (self.gdepth[i]))
            F.write(' %11.3f %11.3f %11.3f %11.3f %11d %11d %11d\n' % (self.LON[0], self.LON[1], self.LAT[0], self.LAT[1], self.nx, self.ny, self.nz))
            for i in range(self.nz+1):
                F.write(' %11.3f\n' % (self.rdepth[i]))
            F.write(' %11d %11d\n' % (self.nmax, self.no_params))
            F.write(' '.join(self.grid.astype(str)))

    def extract_1D(self, lon, lat, values, ref=None, method='linear'):
        """
        Extracts arbitrary 1D vertical seismic velocity profiles from an
        inversion model.

        Input::
            * lon - float or list of desired longitude locations.
            * lat - float or list of desired latitude locations.
            * values - array of inversion values corresponding to the grid
              class.
            * ref [optional] - the reference model to convert velocity
              perturbations to absolute values, if desired. Needs to be a 2D
              array with depth in first column and velocity in second column.
            * method [optional] - method used to interpolate values to location
        """
        if type(lon) != np.ndarray:
            if type(lon) == float:
                lon = np.array([lon])
            else:
                lon = np.array(lon)
        if type(lat) != np.ndarray:
            if type(lat) == float:
                lat = np.array([lat])
            else:
                lat = np.array(lat)

        # Create underlying model meshgrid
        Y = np.arange(self.LAT.max()-self.dy/2., self.LAT.min(), -self.dy)
        X = np.arange(self.LON.min()+self.dx/2., self.LON.max(),  self.dx)
        grid_y, grid_z, grid_x = np.meshgrid(Y, (self.rdepth[:-1] + self.rdepth[1:])/2., X)
        pts = np.vstack((grid_x.reshape(np.prod(grid_x.shape)),
                         grid_y.reshape(np.prod(grid_y.shape)),
                         grid_z.reshape(np.prod(grid_z.shape)))).T
        values = values.reshape(np.prod(values.shape))

        # Create 1D profile meshgrid
        Z = (self.rdepth[:-1] + self.rdepth[1:])/2.
        profiles = np.zeros((lon.shape[0], Z.shape[0]))

        for i in range(lon.shape[0]):
            smp_y, smp_z, smp_x = np.meshgrid(lat[i], Z, lon[i])

            # Sub-sample threshold levels to clustered grid
            profile = griddata(pts, values, (smp_x, smp_y, smp_z), method=method)
            profile = profile.flatten()
            profiles[i] = profile

        if np.any(ref):
            ref = griddata(ref[:,0], ref[:,1], Z, method='linear')
            profiles = (profiles/100. + 1.) * ref

        return profiles, Z

if __name__=='__main__':
    # Needed files: param *locsol.it3 AK135.txt
    phase = 'S'
    lons = [135., 140.21]
    lats = [-25., -15.43]

    # Extract the 1D profiles
    g = grid()
    v = np.loadtxt(phase+'locsol.it3', skiprows=2)
    ak135 = np.loadtxt('AK135.txt', skiprows=3)
    for i in range(ak135.shape[0]):
        if ak135[i,0] == ak135[i-1,0]:
            ak135[i,0] += 0.001
    if phase == 'P':
        ak135 = ak135[:,[True, False, True, False, False, False, False]]
    else:
        ak135 = ak135[:,[True, False, False, True, False, False, False]]
    profiles, Z = g.extract_1D(lons, lats, v, ref=ak135)

    # Demonstrate recovery
    from matplotlib import pyplot as plt
    labels = ['(%.2f$^\circ$E, %.2f$^\circ$S)' % (lons[i], -lats[i]) for i in range(len(lons))]
    for i in range(len(lons)):
        plt.plot(profiles[i], Z, label=labels[i])
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
