#!/bin/env python
"""
Description:
    Implements migration algorithm as described in Frassetto et al. (2010):
      Improved imaging with phase-weighted common conversion point stacks
      of receiver functions (GJI)

References:

CreationDate:   3/15/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     3/15/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
import logging
import click
from collections import defaultdict

import numpy as np

import pkg_resources
import matplotlib as mpl
import matplotlib.pyplot as plt
from obspy import read_inventory, read_events, UTCDateTime as UTC
from obspy.clients.fdsn import Client
from rf import read_rf, RFStream
from rf import get_profile_boxes, iter_event_data, IterMultipleComponents
from rf.util import _add_processing_info, direct_geodetic

from scipy.spatial import cKDTree
from scipy.interpolate import interp1d
from scipy.signal import hilbert
from mpi4py import MPI

# pylint: disable=invalid-name

logging.basicConfig()
log = logging.getLogger('migration')


# define utility functions
def rtp2xyz(r, theta, phi):
    """Convert spherical to cartesian coordinates

    :param r: [description]
    :type r: [type]
    :param theta: [description]
    :type theta: [type]
    :param phi: [description]
    :type phi: [type]
    :return: [description]
    :rtype: [type]
    """
    xout = np.zeros((r.shape[0], 3))
    rst = r * np.sin(theta)
    xout[:, 0] = rst * np.cos(phi)
    xout[:, 1] = rst * np.sin(phi)
    xout[:, 2] = r * np.cos(theta)
    return xout
# end func

def xyz2rtp(x, y, z):
    """Convert cartesian to spherical coordinates

    :param x: [description]
    :type x: [type]
    :param y: [description]
    :type y: [type]
    :param z: [description]
    :type z: [type]
    :return: [description]
    :rtype: [type]
    """
    rout = np.zeros((x.shape[0], 3))
    tmp1 = x * x + y * y
    tmp2 = tmp1 + z * z
    rout[0] = np.sqrt(tmp2)
    rout[1] = np.arctan2(np.sqrt(tmp1), z)
    rout[2] = np.arctan2(y, x)
    return rout
# end func

class Geometry:
    def __init__(self, start_lat_lon, azimuth, lengthkm, nx, widthkm, ny, depthkm, nz, debug=False):
        self._start_lat_lon = np.array(start_lat_lon)
        assert self._start_lat_lon.shape == (2,), 'start lat-lon should be a list of length 2'

        self._azimuth = azimuth
        self._length = lengthkm
        self._width = widthkm
        self._depth = depthkm
        self._nx = nx
        self._ny = ny
        self._nz = nz

        self._ortho = (self._azimuth + 90) % 360 # orthogonal to azimuth
        self._earth_radius = 6371 #km
        self._debug = debug
        # Generate four sets of grids:
        # 1. Lon-Lat-Depth grid, with slowest to fastest index in that order
        # 2. Cartesian axis-aligned x-y-z grid, with slowest to fastest index in
        #    that order, starting from start_lat_lon. Note that this is a local
        #    coordinate system that does not account for spherical curvature and
        #    used for plotting purposes alone
        # 3. Spherical grid in Cartesian coordinates that accounts for spherical
        #    curvature and is used internally for all nearest neighbour calculations
        # 4. Cell-centre coordinates for 3
        self._glon, self._glat, self._gz, self._gxaa, self._gyaa, self._gzaa, \
            self._gxs, self._gys, self._gzs, self._gxsc, self._gysc, self._gzsc = self.generateGrids()

        # Compute centres of depth-node pairs
        self._gzaac = (self._gzaa[0,0,1:] + self._gzaa[0,0,:-1])/2.
    # end func

    def generateGrids(self):
        # Start mesh generation==============================================
        sll = self._start_lat_lon

        result = []
        resultCart = []
        dx = self._length / float(self._nx - 1)
        dy = self._width / float(self._ny - 1)
        dz = self._depth / float(self._nz - 1)

        runLengthll = sll
        runWidthll = sll
        # cx = cy = cz = 0
        for ix in range(self._nx):
            runWidthll = runLengthll
            for iy in range(self._ny):
                for iz in range(self._nz):
                    result.append([runWidthll[1], runWidthll[0], iz * dz])
                    resultCart.append([ix * dx, iy * dy, iz * dz])
                # end for
                runWidthll = direct_geodetic(runLengthll, self._ortho, iy * dy)
            # end for
            runLengthll = direct_geodetic(runLengthll, self._azimuth, dx)
        # end for
        result = np.array(result).reshape(self._nx, self._ny, self._nz, 3)
        resultCart = np.array(resultCart).reshape(self._nx, self._ny, self._nz, 3)

        glon = result[:, :, :, 0].reshape(self._nx, self._ny, self._nz)
        glat = result[:, :, :, 1].reshape(self._nx, self._ny, self._nz)

        # Create local cartesian axis-aligned grids
        # Naming convention (Grid [XYZ] Axis-Aligned)
        gxaa = resultCart[:, :, :, 0].reshape(self._nx, self._ny, self._nz)
        gyaa = resultCart[:, :, :, 1].reshape(self._nx, self._ny, self._nz)
        gzaa = resultCart[:, :, :, 2].reshape(self._nx, self._ny, self._nz)

        # Create cartesian mesh with spherical curvature
        # Naming convention (Grid [XYZ] Spherical)
        ts = (90 - glat.flatten()) / 180. * np.pi
        ps = glon.flatten() / 180. * np.pi
        rs = (self._earth_radius - gzaa.flatten()) * np.ones(ts.shape)
        rtps = np.array([rs, ts, ps]).T
        xyzs = rtp2xyz(rtps[:, 0], rtps[:, 1], rtps[:, 2])
        xyzs = np.array(xyzs).reshape(self._nx, self._ny, self._nz, 3)
        gxs = xyzs[:, :, :, 0].reshape(self._nx, self._ny, self._nz)
        gys = xyzs[:, :, :, 1].reshape(self._nx, self._ny, self._nz)
        gzs = xyzs[:, :, :, 2].reshape(self._nx, self._ny, self._nz)

        if(self._debug):
            print(np.min(gxs.flatten()), np.max(gxs.flatten()))
            print(np.min(gys.flatten()), np.max(gys.flatten()))
            print(np.min(gzs.flatten()), np.max(gzs.flatten()))

        # Compute cell-centre coordinates
        # Naming convention (Grid [XYZ] Spherical Centre)
        gxsc = (gxs[:-1, :-1, :-1] + gxs[1:, 1:, 1:]) / 2.
        gysc = (gys[:-1, :-1, :-1] + gys[1:, 1:, 1:]) / 2.
        gzsc = (gzs[:-1, :-1, :-1] + gzs[1:, 1:, 1:]) / 2.

        if(self._debug):
            print('\n')
            print(np.min(gxsc.flatten()), np.max(gxsc.flatten()))
            print(np.min(gysc.flatten()), np.max(gysc.flatten()))
            print(np.min(gzsc.flatten()), np.max(gzsc.flatten()))

        return glon, glat, gzaa, gxaa, gyaa, gzaa, gxs, gys, gzs, gxsc, gysc, gzsc
    # end func
# end class

class Migrate:
    def __init__(self, geometry, stream, debug=False, output_folder='/tmp'):
        assert isinstance(geometry, Geometry), 'Must be an instance of class Geometry..'
        self._geometry = geometry

        assert isinstance(stream, RFStream), 'Must be an instance of class RFStream..'
        self._stream = stream

        self._debug = debug
        self._output_folder = output_folder
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Initialize MPI
        self._comm = MPI.COMM_WORLD
        self._nproc = self._comm.Get_size()
        self._chunk_index = self._comm.Get_rank()

        self._ppDict = defaultdict(list) # dictionary for piercing point results
        self._proc_izs = defaultdict(list) # depth indices that each process works on
        self._treeDict = {} # dictionary for Kd-trees for each depth layer
        self._d2tIO = None
        self._iphaseTraceDataList = []

        # Create depth-to-time interpolation object
        # rp = 'rf'
        # rpf = '/'.join(('data', 'ak135.dat')) # find where ak135.dat is located
        #fp = pkg_resources.resource_stream(rp, rpf)
        #fn = fp.name
        #fp.close()

        velocity_model_file = os.path.join(os.path.split(__file__)[0], 'models', 'iasp91.dat')

        m = np.loadtxt(velocity_model_file)
        dlim = m[:, 0] < 2800 # don't need data past 2800 km

        depths = m[:, 0][dlim] # depths in km
        s = 1. / m[:, 2][dlim] # slowness in s/km
        # Integrate slowness with respect to distance to get travel times.
        # TODO: should we use Dix's interval velocity?
        times = np.cumsum(np.diff(depths) * (s[0:-1] + s[1:]) / 2.)
        times = np.insert(times, 0, 0)
        self._d2tIO = interp1d(depths, times)

        # split workload
        self.__split_work()

        # compute instantaneous phase
        for t in self._stream:
            analytic = hilbert(t.data)
            angle = np.angle(analytic)
            iPhase = np.exp(1j * angle)
            self._iphaseTraceDataList.append(iPhase)
        #end for

        log.info('Computed instantaneous phases.')
    # end func

    def __split_work(self):
        """
        Splits up workload over n processors
        """

        if (self._chunk_index == 0):
            count = 0
            for iproc in np.arange(self._nproc):
                for iz in np.arange(np.divide(self._geometry._gzaac.shape[0], self._nproc)):
                    self._proc_izs[iproc].append(count)
                    count += 1
            # end for

            for iproc in np.arange(np.mod(self._geometry._gzaac.shape[0], self._nproc)):
                self._proc_izs[iproc].append(count)
                count += 1
        # end if

        # broadcast workload to all procs
        self._proc_izs = self._comm.bcast(self._proc_izs, root=0)
        if (self._chunk_index == 0): log.info(' Distributing workload over %d processors..' % (self._nproc))

        if(self._debug):
            print('proc: %d, %d depth values\n========='%(self._chunk_index,
                                                   len(self._proc_izs[self._chunk_index])))
            for iz in self._proc_izs[self._chunk_index]:
                print(iz, self._geometry._gzaac[iz])
        # end if
    # end func

    def __generatePiercingPoints(self):
        for iz in self._proc_izs[self._chunk_index]:
            z = self._geometry._gzaac[iz]
            ppoints = self._stream.ppoints(z)
            self._ppDict[iz] = ppoints
            #print iz, z, len(ppoints)
        # end for

        # Gather all results on proc 0
        ppDictList = self._comm.gather(self._ppDict, root=0)

        if(self._chunk_index==0):
            self._ppDict = defaultdict(list)
            for ip in np.arange(self._nproc):
                for k in ppDictList[ip].keys():
                    self._ppDict[k] = ppDictList[ip][k]

                    #print len(ppDictList[ip][k])
                # end for
            # end for

            if(self._debug):
                f = open(os.path.join(self._output_folder, 'pp_parallel.txt'), 'w+')
                for k in sorted(self._ppDict.keys()):
                    f.write('%f\n'%(self._geometry._gzaac[k]))
                    pp = self._ppDict[k]
                    for i in pp:
                        f.write('%f %f\n'%(i[0], i[1]))
                f.close()
            #end if

            # Create xyz nodes for each depth value
            for k in self._ppDict.keys():
                ts = (90 - self._ppDict[k][:, 0]) / 180. * np.pi
                ps = self._ppDict[k][:, 1] / 180. * np.pi

                z = self._geometry._gzaac[k]
                rs = (self._geometry._earth_radius - z) * np.ones(ts.shape)
                rtps = np.array([rs, ts, ps]).T

                xyzs = rtp2xyz(rtps[:, 0], rtps[:, 1], rtps[:, 2])
                self._treeDict[k] = xyzs
            # end for
        # end if

        if(self._chunk_index==0): log.info(' Broadcasting Kd-tree nodes..')
        # broadcast xyz nodes to all procs
        self._treeDict = self._comm.bcast(self._treeDict, root=0)
        
        # Create Kd-trees for each depth value
        for k in sorted(self._treeDict.keys()):
            self._treeDict[k] = cKDTree(self._treeDict[k])
        
        if(self._debug):
            f = open(os.path.join(self._output_folder, 'kd_parallel.%02d.txt'%(self._chunk_index)), 'w+')
            for k in sorted(self._treeDict.keys()):
                f.write('%f\n'%(self._geometry._gzaac[k]))
                for i in self._treeDict[k].data:
                    f.write('%f %f %f\n' % (i[0], i[1], i[2]))
            f.close()
        # end if
    # end func

    def execute(self):
        if(self._chunk_index==0): log.info(' Generating Piercing Points..')
        self.__generatePiercingPoints()

        if (self._chunk_index == 0): log.info(' Stacking amplitudes..')

        vol = np.zeros(self._geometry._gxsc.shape)
        cz = np.zeros(self._geometry._gxsc.shape, dtype=np.complex_)
        volHits = np.zeros(self._geometry._gxsc.shape)
        
        times = self._stream[0].times() - 25
        for ix in range(self._geometry._nx - 1):
            for iy in range(self._geometry._ny - 1):
                for iz in self._proc_izs[self._chunk_index]:
                    z = self._geometry._gzaac[iz]
                    t = self._treeDict[iz]

                    ids = t.query_ball_point([self._geometry._gxsc[ix, iy, iz],
                                              self._geometry._gysc[ix, iy, iz],
                                              self._geometry._gzsc[ix, iy, iz]], r=20)

                    if (len(ids) == 0):
                        continue
                    # end if

                    ct = self._d2tIO(z)
                    tidx = np.argmin(np.fabs(times - ct))
                    for i, si in enumerate(ids):
                        # print tidx*(1./stream[si].stats.sampling_rate)-25, d2tIO(z)

                        vol[ix, iy, iz] += self._stream[si].data[tidx]
                        cz[ix, iy, iz] += self._iphaseTraceDataList[si][tidx]

                        volHits[ix, iy, iz] += 1.
                    # end for

                    if(volHits[ix, iy, iz]>0):
                        cz[ix, iy, iz] /= volHits[ix, iy, iz]
                    
                    cz[ix, iy, iz] = np.abs(cz[ix, iy, iz])
                # end for
            # end for
            #print ix
        # end for
        
        cz = cz.astype(np.float64)

        # Sum all on master proc
        if(self._chunk_index==0):
            totalCz = np.zeros(self._geometry._gxsc.shape)
            totalVol = np.zeros(self._geometry._gxsc.shape)
            totalVolHits = np.zeros(self._geometry._gxsc.shape)
        else:
            totalCz = None
            totalVol = None
            totalVolHits = None
        # end if

        self._comm.Reduce([cz, MPI.DOUBLE], [totalCz, MPI.DOUBLE],
                          op=MPI.SUM, root=0)
        self._comm.Reduce([vol, MPI.DOUBLE], [totalVol, MPI.DOUBLE],
                          op=MPI.SUM, root=0)
        self._comm.Reduce([volHits, MPI.DOUBLE], [totalVolHits, MPI.DOUBLE],
                          op=MPI.SUM, root=0)

        if(self._chunk_index==0):
            np.savetxt(os.path.join(self._output_folder, 'cz.txt'), totalCz.flatten())
            np.savetxt(os.path.join(self._output_folder, 'vol.txt'), totalVol.flatten())
            np.savetxt(os.path.join(self._output_folder, 'volHits.txt'), totalVolHits.flatten())
            np.savetxt(os.path.join(self._output_folder, 'gxaa.txt'), self._geometry._gxaa.flatten())
            np.savetxt(os.path.join(self._output_folder, 'gyaa.txt'), self._geometry._gyaa.flatten())
            np.savetxt(os.path.join(self._output_folder, 'gzaa.txt'), self._geometry._gzaa.flatten())
            np.savetxt(os.path.join(self._output_folder, 'glon.txt'), self._geometry._glon.flatten())
            np.savetxt(os.path.join(self._output_folder, 'glat.txt'), self._geometry._glat.flatten())
    # end func
# end class


@click.command()
@click.option('--start-lat-lon', type=(float, float), required=True, help='Start latitude, longitude (degrees)')
@click.option('--azimuth', type=float, required=True, help='document me')
@click.option('--dimensions', type=(float, float, float), required=True,
              help='(length, width, depth) of the volume, in km')
@click.option('--num-cells', type=(int, int, int), required=True,
              help='Number of discrete cells in each dimension of the volume')
@click.option('--debug', type=bool, is_flag=True, default=False)
@click.argument('rf-h5-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-folder', type=click.Path(file_okay=False), required=True)
def main(rf_h5_file, output_folder, start_lat_lon, azimuth, dimensions, num_cells, debug):
    """Perform 3D migration of RFs to volumetric space, stacking RF amplitudes in each cell.

    Example usage:
        python rf_3dmigrate.py --start-lat-lon -17.4 132.9 --azimuth 80 --dimensions 1000 450 75 \
            --num-cells 100 45 375 /g/data/ha3/am7399/shared/OA-ZRT-R-cleaned.h5 /g/data/ha3/am7399/shared/OA_piercing

    :param rf_h5_file: Source file containing receiver functions
    :type rf_h5_file: str or Path
    :param output_folder: Folder in which to output results
    :type output_folder: str or Path
    """

    s = read_rf(rf_h5_file, 'H5')

    g = Geometry(start_lat_lon=start_lat_lon, azimuth=azimuth,
                 lengthkm=dimensions[0], nx=num_cells[0],
                 widthkm=dimensions[1], ny=num_cells[1],
                 depthkm=dimensions[2], nz=num_cells[2], debug=debug)

    m = Migrate(geometry=g, stream=s, debug=False, output_folder=output_folder)
    m.execute()
# end


if __name__ == "__main__":
    # call main function
    main()  # pylint: disable=no-value-for-parameter
