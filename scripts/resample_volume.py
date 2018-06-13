#!/bin/env python
"""
Description:
    The Resample class reads in a xyz file, containing longitude, latitude, depth and value
    quadruplets and resamples these values on a finer grid based on inverse distance weighted
    interpolation (IDW)
   
References: https://en.wikipedia.org/wiki/Inverse_distance_weighting
 
CreationDate:   6/5/18
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     6/5/18   RH

"""

import os
import numpy as np
from scipy.spatial import cKDTree
from scipy.interpolate import interp1d
from pyevtk.hl import gridToVTK
from multiprocessing import Pool

# define utility functions
def rtp2xyz(r, theta, phi):
    xout = np.zeros((r.shape[0], 3))
    rst = r * np.sin(theta);
    xout[:, 0] = rst * np.cos(phi)
    xout[:, 1] = rst * np.sin(phi)
    xout[:, 2] = r * np.cos(theta)
    return xout
# end func

def xyz2rtp(x, y, z):
    rout = np.zeros((x.shape[0], 3))
    tmp1 = x * x + y * y
    tmp2 = tmp1 + z * z
    rout[:, 0] = np.sqrt(tmp2)
    rout[:, 1] = np.arctan2(np.sqrt(tmp1), z)
    rout[:, 2] = np.arctan2(y, x)
    return rout
# end func

class Resample:
    def __init__(self, fname, lonf, latf, zf, nn=1, p=4):
        """
        :param fname: xyz filename
        :param lonf: factor by which to increase longitudinal sampling (>1)
        :param latf: factor by which to increase latitudinal sampling (>1)
        :param zf: factor by which to increase z sampling (>1)

        """
        if not os.path.isfile(fname):
            raise RuntimeError('File %s not found..'%fname)
        if(lonf < 1 or latf < 1 or zf < 1):
            raise RuntimeError('Invalid sampling factor..')

        self.RADIUS = 6371
        self._fname = fname
        self._lonf = lonf
        self._latf = latf
        self._zf = zf

        # file format [lon, lat, z(depth km), perturbation]
        data = np.loadtxt(self._fname)
        self._gridLLZ = data[:, :3]
        self._values = data[:, -1]

        self._gridLLZ[:, 2] = self.RADIUS - self._gridLLZ[:, 2]

        # get mesh dimensions
        indices = np.where(self._gridLLZ[0, 0] == self._gridLLZ[:, 0])
        self._nx = indices[0][1] - indices[0][0]
        indices = np.where(self._gridLLZ[0, 1] == self._gridLLZ[::self._nx, 1])
        self._ny = indices[0][1] - indices[0][0]
        self._nz = np.shape(self._gridLLZ[::self._nx * self._ny, 2])[0]

        # convert to Cartesian coordinates
        self._gridXYZ  = rtp2xyz(self._gridLLZ[:, 2],
                                 np.radians(90-self._gridLLZ[:, 1]),
                                 np.radians(self._gridLLZ[:, 0]))

        # create kdTree
        self._gridTree = cKDTree(self._gridXYZ)

        # create denser mesh
        self._rnx = int(self._nx * self._lonf)
        self._rny = int(self._ny * self._latf)
        self._rnz = int(self._nz * self._zf)

        lono = self._gridLLZ[:self._nx, 0] # original lons
        lato = self._gridLLZ[:self._nx*self._ny:self._nx, 1] # original lats
        zo = self._gridLLZ[::self._nx*self._ny, 2] # original zs

        # preserve distribution of nodes, especially for z nodes,
        # which become sparser as we go deeper
        lonio = interp1d(np.linspace(0, len(lono), len(lono)), lono)
        latio = interp1d(np.linspace(0, len(lato), len(lato)), lato)
        zio = interp1d(np.linspace(0, len(zo), len(zo)), zo)

        lonr = lonio(np.linspace(0, len(lono), len(lono) * self._lonf))
        latr = latio(np.linspace(0, len(lato), len(lato) * self._latf))
        zr = zio(np.linspace(0, len(zo), len(zo) * self._zf))

        self._rgridLon, self._rgridLat, self._rgridZ = np.meshgrid(lonr, latr, zr, indexing='ij')

        # transpose axes to conform to original ordering
        self._rgridLon = np.transpose(self._rgridLon, (2, 1, 0))
        self._rgridLat = np.transpose(self._rgridLat, (2, 1, 0))
        self._rgridZ = np.transpose(self._rgridZ, (2, 1, 0))

        self._rgridXYZ = rtp2xyz(self._rgridZ.flatten(),
                                 np.radians(90 - self._rgridLat.flatten()),
                                 np.radians(self._rgridLon.flatten()))

        # Resample
        # query Kd-tree instance to retrieve distances and
        # indices of k nearest neighbours

        img = []
        vals = self._values.flatten()

        rg = np.reshape(self._rgridXYZ, (self._rnx, self._rny, self._rnz, 3))

        for iz in range(self._rnz):

            xyzs = rg[:, :, iz, :]
            xyzs = xyzs.reshape(-1, xyzs.shape[-1])

            d, l = self._gridTree.query(xyzs, k=nn)

            if (nn == 1):
                # extract nearest neighbour values
                img.append(self._values.flatten()[l])
            else:
                result = np.zeros((self._rnx * self._rny))
                # field values are directly assigned for coincident locations
                coincidentValIndices = d[:, 0] == 0
                result[coincidentValIndices] = vals[l[coincidentValIndices, 0]]

                # perform idw interpolation for non-coincident locations
                idwIndices = d[:, 0] != 0
                w = np.zeros(d.shape)
                w[idwIndices, :] = 1. / np.power(d[idwIndices, :], p)

                result[idwIndices] = np.sum(w[idwIndices, :] * vals[l[idwIndices, :]], axis=1) / \
                                     np.sum(w[idwIndices, :], axis=1)

                img.append(result)
            # end if
        # end for

        self._img = np.array(img).transpose(1, 0).flatten()
    # end func

    def GetOriginalResolution(self):
        """

        :return: nx, ny, nz
        """
        return self._nx, self._ny, self._nz
    # end func

    def WriteVTK(self, fname):
        """

        :return:
        """
        lons, lats, ds, xs, ys, zs =  self._rgridLon.flatten(), self._rgridLat.flatten(), \
                                      self.RADIUS - self._rgridZ.flatten(), np.array(self._rgridXYZ[:, 0]), \
                                      np.array(self._rgridXYZ[:, 1]), np.array(self._rgridXYZ[:, 2])

        cv = np.reshape(self._img, (self._rnx, self._rny, self._rnz), order='F')
        cd = np.reshape(ds, (self._rnx, self._rny, self._rnz), order='F')
        cv = (cv[1:, 1:, 1:] + cv[:-1, :-1, :-1]) / 2.
        cd = (cd[1:, 1:, 1:] + cd[:-1, :-1, :-1]) / 2.

        cv = np.reshape(cv, cv.shape, order='F')
        cd = np.reshape(cd, cd.shape, order='F')

        gridToVTK(fname,
                  np.reshape(xs, (self._rnx, self._rny, self._rnz), order='F'),
                  np.reshape(ys, (self._rnx, self._rny, self._rnz), order='F'),
                  np.reshape(zs, (self._rnx, self._rny, self._rnz), order='F'),
                  cellData={'Amplitude': cv, 'd': cd})

        # np.savetxt('/tmp/tomo.txt', np.array([lons, lats, zs, values]).T)
        # end func
# end class

def main():
    """
    define main function
    :return:
    """

    flon = 5
    flat = 5
    fz = 5

    rs = Resample('/home/rakib/work/pst/tomo/P-wave_1x1.it3.xyz', flon, flat, fz, nn=50, p=0.5)
    rs.WriteVTK('/tmp/tomoRadVariableSmoothing')

    return
# end func


# =============================================
# Quick test
# =============================================
if __name__ == "__main__":
    # call main function
    main()