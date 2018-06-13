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
from netCDF4 import Dataset
import click
import os.path as path

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

    def WriteVTK(self, fnamestem):
        """

        :return:
        """
        lons, lats, ds, xs, ys, zs =  self._rgridLon.flatten(), self._rgridLat.flatten(), \
                                      self.RADIUS - self._rgridZ.flatten(), np.array(self._rgridXYZ[:, 0]), \
                                      np.array(self._rgridXYZ[:, 1]), np.array(self._rgridXYZ[:, 2])

        v = np.reshape(self._img, (self._rnx, self._rny, self._rnz), order='F')
        d = np.reshape(ds, (self._rnx, self._rny, self._rnz), order='F')
        cv = (v[1:, 1:, 1:] + v[:-1, :-1, :-1]) / 2.
        cd = (d[1:, 1:, 1:] + d[:-1, :-1, :-1]) / 2.

        cv = np.reshape(cv, cv.shape, order='F')
        cd = np.reshape(cd, cd.shape, order='F')

        gridToVTK(fnamestem,
                  np.reshape(xs, (self._rnx, self._rny, self._rnz), order='F'),
                  np.reshape(ys, (self._rnx, self._rny, self._rnz), order='F'),
                  np.reshape(zs, (self._rnx, self._rny, self._rnz), order='F'),
                  cellData={'p': cv, 'Depth': cd})
    # end func

    def WriteNetcdf(self, fnamestem):
        """

        :return:
        """
        lons, lats, ds, xs, ys, zs =  self._rgridLon.flatten(), self._rgridLat.flatten(), \
                                      self.RADIUS - self._rgridZ.flatten(), np.array(self._rgridXYZ[:, 0]), \
                                      np.array(self._rgridXYZ[:, 1]), np.array(self._rgridXYZ[:, 2])

        lons = np.reshape(lons, (self._rnx, self._rny, self._rnz), order='F')
        lats = np.reshape(lats, (self._rnx, self._rny, self._rnz), order='F')
        d = np.reshape(ds, (self._rnx, self._rny, self._rnz), order='F')
        v = np.reshape(self._img, (self._rnx, self._rny, self._rnz), order='F')

        fn = fnamestem + '.nc'

        root_grp = Dataset(fn, 'w', format='NETCDF4')
        root_grp.description = path.basename(fnamestem)

        # Dimensions
        root_grp.createDimension('lon', self._rnx)
        root_grp.createDimension('lat', self._rny)
        root_grp.createDimension('depth', self._rnz)

        # Variables
        lon = root_grp.createVariable('lon', 'f4', ('lon',))
        lat = root_grp.createVariable('lat', 'f4', ('lat',))
        depth = root_grp.createVariable('depth', 'f4', ('depth',))
        amplitude = root_grp.createVariable('p', 'f4', ('lon', 'lat', 'depth'))

        lon[:] = lons[:, 0, 0]
        lat[:] = lats[0, :, 0]
        depth[:] = d[0, 0, :]
        amplitude[:, :, :] = v

        root_grp.close()
    # end func

    def WriteSGrid(self, fnamestem):
        def write_header(hfn, dfn, mname):
            hlines = [  'GOCAD SGrid 1',
                        'HEADER {',
                        'name:%s'%(mname),
                        '*volume:true',
                        'ascii:on',
                        '*painted:on',
                        '*painted*variable:p',
                        'last_selected_folder:Volumes',
                        '*volume*grid:true',
                        '*volume*solid:true',
                        '*cube*solid:false',
                        '*cube*grid:false',
                        '*volume*points:false',
                        '}',
                        'AXIS_N %d %d %d' % (self._rnx, self._rny, self._rnz),
                        'PROP_ALIGNMENT POINTS',
                        'ASCII_DATA_FILE %s' % path.basename(dfn),
                        'PROPERTY 1 "p"',
                        'PROPERTY_CLASS 1 "p"',
                        'END]' ]

            fo = open(hfn, 'w+')
            for l in hlines: fo.write(l+'\n')
            fo.close()
        # end func

        def write_data(dfn):
            lons, lats, ds = self._rgridLon.flatten(), self._rgridLat.flatten(), \
                             self.RADIUS - self._rgridZ.flatten()

            kxyz, jxyz, ixyz = np.meshgrid(np.arange(self._rnz), np.arange(self._rny),
                                           np.arange(self._rnx), indexing='ij')

            np.savetxt(dfn, np.array([lons, lats, ds*1e3, self._img,
                                      ixyz.flatten(), jxyz.flatten(),
                                      kxyz.flatten()]).T, fmt='%4.4f %4.4f %4.4f %5.7f %d %d %d')
        # end func

        hfn = fnamestem + '.sg'
        dfn = fnamestem + '.ascii.txt'
        modelname = path.basename(fnamestem)

        write_header(hfn, dfn, modelname)
        write_data(dfn)
    # end func
# end class

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input', required=True,
                type=click.Path(exists=True))
@click.option('--output-file-stem', default='none', type=str,
              help="Stem of output file; the default, 'none', uses the input file name without the extension")
@click.option('--output-type', default='sgrid', type=click.Choice(['sgrid', 'vtk', 'netcdf']),
              help="Output file-type; default 'sgrid'")
@click.option('--resampling-factors', default='5x5x5', type=str,
              help="Resampling factors for the longitudinal, latitudinal and depth axes. The default value of "
                   "'5x5x5' implies the output mesh will be 5 times denser in the longitudinal, latitudinal and "
                   "depth directions, respectively")
@click.option('--nearest-neighbours', default=50, type=int,
              help="Number of nearest neighbours to include in Inverse-Distance-Weighted (IDW) interpolation.")
@click.option('--power-parameter', default=0.5, type=float,
              help="Power parameter, which determines the relative influence of near and far neighbours during "
                   "interpolation.")
def process(input, output_file_stem, output_type, resampling_factors,
            nearest_neighbours, power_parameter):
    """
    Script for super- or sub-sampling a 3D volume using inverse-distance-weighted interpolation.

    INPUT: Path to xyz file, with nodal values on a regular 3D grid\n

    Example Usage: python resample_volume.py P-wave_1x1.it3.xyz --resampling-factors 2x2x2 --output-file-stem /tmp/tomo
    """
    try:
        flon, flat, fz = np.float_(resampling_factors.lower().split('x'))
    except:
        raise RuntimeError('Invalid resampling factors')
    # end try

    if(output_file_stem=='none'): output_file_stem = path.splitext(input)[0]

    rs = Resample(input, flon, flat, fz, nn=nearest_neighbours, p=power_parameter)

    if(output_type=='sgrid'):
        rs.WriteSGrid(output_file_stem)
    elif(output_type=='vtk'):
        rs.WriteVTK(output_file_stem)
    else:
        rs.WriteNetcdf(output_file_stem)
    # end if
# end func

if __name__ == "__main__":
    process()