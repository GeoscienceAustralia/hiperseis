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
from seismic.receiver_fn import rf_util, rf_corrections
import obspy.geodetics.base as geodetics
import numpy as np
from rf import read_rf, RFStream
from rf.util import direct_geodetic, DEG2KM
from obspy.taup import TauPyModel
import h5py
from tqdm import tqdm

from scipy.spatial import cKDTree
from scipy.interpolate import splev, splrep, sproot, interp1d, InterpolatedUnivariateSpline
from scipy import optimize
from scipy.integrate import simps as simpson
from scipy.signal import hilbert
from mpi4py import MPI


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
        self._glon, self._glat, self._gz, self._gxl, self._gyl, self._gzl, \
            self._gxs, self._gys, self._gzs, self._gxsc, self._gysc, self._gzsc = self.generateGrids()

        # Compute centres of depth-node pairs
        self._gzlc = (self._gzl[0,0,1:] + self._gzl[0,0,:-1])/2.
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
        gxl = resultCart[:, :, :, 0].reshape(self._nx, self._ny, self._nz)
        gyl = resultCart[:, :, :, 1].reshape(self._nx, self._ny, self._nz)
        gzl = resultCart[:, :, :, 2].reshape(self._nx, self._ny, self._nz)

        # Create cartesian mesh with spherical curvature
        # Naming convention (Grid [XYZ] Spherical)
        ts = (90 - glat.flatten()) / 180. * np.pi
        ps = glon.flatten() / 180. * np.pi
        rs = (self._earth_radius - gzl.flatten()) * np.ones(ts.shape)
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

        return glon, glat, gzl, gxl, gyl, gzl, gxs, gys, gzs, gxsc, gysc, gzsc
    # end func
# end class

class Migrate:
    def __init__(self, rf_filename, debug=False, output_file=None, earth_radius=6371):
        self._rf_filename = rf_filename
        self._debug = debug
        self._output_file = output_file
        self._earth_radius = earth_radius #km

        # Initialize MPI
        self._comm = MPI.COMM_WORLD
        self._nproc = self._comm.Get_size()
        self._rank = self._comm.Get_rank()

        self._p_traces = None # primary RF traces
        self._iphase_traces = None # instantaneous phase traces
        self._stations = defaultdict(list)
        self._process_streams()

        log.info('Computed instantaneous phases.')
    # end func

    def _process_streams(self, primary_component='R', model='iasp91', phase='P'):
        proc_hkeys = None
        if(self._rank == 0):
            hkeys = rf_util.get_hdf_keys(self._rf_filename)
            assert len(hkeys), 'No hdf5-groups found in file {}. Aborting..'.format(self._rf_filename)

            # split workload
            proc_hkeys = rf_util.split_list(hkeys, self._nproc)
        # end if
        # broadcast workload to all procs
        proc_hkeys = self._comm.bcast(proc_hkeys, root=0)

        # read local traces

        self._p_traces = RFStream()
        #pbar = tqdm(total=len(proc_hkeys[self._rank]))
        for hkey in proc_hkeys[self._rank]:
            #pbar.update()
            #pbar.set_description('rank {}: loading {} '.format(self._rank, hkey))

            traces = read_rf(self._rf_filename, format='h5', group='waveforms/%s'%(hkey))

            # select primary component only
            p_traces = traces.select(component=primary_component)
            assert len(p_traces), 'No {} component found in RF stream for {}. Aborting..'.format(primary_component,
                                                                                                 hkey)
            min_slope_ratio = 5
            if (min_slope_ratio > 0):
                p_traces = RFStream([tr for tr in p_traces \
                                     if tr.stats.slope_ratio > min_slope_ratio])
            # end if

            if(len(p_traces) == 0): continue

            has_reverberations = rf_corrections.has_reverberations(p_traces)
            if(has_reverberations):
                print(hkey)
                p_traces = rf_corrections.apply_reverberation_filter(p_traces)
            # end if

            p_traces.filter(type='highpass', freq=0.1,
                            corners=1, zerophase=True)

            self._p_traces += p_traces
            self._stations[hkey] = [p_traces[0].stats.station_longitude,
                                    p_traces[0].stats.station_latitude,
                                    p_traces[0].stats.station_elevation,
                                    has_reverberations]
        # end for
        #pbar.close()

        if(not len(self._p_traces)):
            return
        # end if

        # compute instantaneous phases
        self._iphase_traces = self._p_traces.copy()
        for t in self._iphase_traces:
            analytic = hilbert(t.data)
            angle = np.angle(analytic)
            iphase = np.exp(1j * angle)
            t.data = iphase
        #end for

        # compute ray-paths
        taup_model = TauPyModel(model=model)
        vmodel = taup_model.model.s_mod.v_mod

        #pbar = tqdm(total=len(self._p_traces))
        point_data_array = None
        for it, (t, ipt) in enumerate(zip(self._p_traces, self._iphase_traces)):
            #pbar.update()
            desc = 'rank {}: computing p-points for {} '.format(self._rank,
                                                             '.'.join([t.stats.network,
                                                             t.stats.station,
                                                             t.stats.location]))
            #pbar.set_description(desc)
            arrivals = taup_model.get_ray_paths_geo(source_depth_in_km=t.stats.event_depth,
                                                    source_latitude_in_deg=t.stats.event_latitude,
                                                    source_longitude_in_deg=t.stats.event_longitude,
                                                    receiver_latitude_in_deg=t.stats.station_latitude,
                                                    receiver_longitude_in_deg=t.stats.station_longitude,
                                                    phase_list=[phase], resample=True)
            arr = arrivals[0]

            times  = (arr.path['time'][-1] - arr.path['time'])[::-1]
            depths = arr.path['depth'][::-1]
            lats   = arr.path['lat'][::-1]
            lons   = arr.path['lon'][::-1]

            rs     = self._earth_radius - depths
            thetas = np.radians(90 - lats)
            phis   = np.radians(lons)

            xyzs   = rtp2xyz(rs, thetas, phis)

            t2dio  = interp1d(times, depths)
            t2xio  = interp1d(times, xyzs[:,0])
            t2yio  = interp1d(times, xyzs[:,1])
            t2zio  = interp1d(times, xyzs[:,2])

            tck = splrep(times, depths, k=3, s=0)

            NZ = 600
            ZMAX = 150
            dnew = np.linspace(1e-5, ZMAX, NZ)
            tnew = np.zeros(dnew.shape)
            for idx in np.arange(dnew.shape[0]):
                tck_mod = (tck[0], tck[1] - dnew[idx], tck[2])
                tnew[idx] = sproot(tck_mod)[0]
            #end for

            d2tio = interp1d(np.concatenate(([0], dnew)), np.concatenate(([0], tnew)))
            vp = vmodel.evaluate_above(dnew, 'p')
            vs = vmodel.evaluate_above(dnew, 's')

            tps = np.zeros(dnew.shape)
            p = t.stats.slowness / DEG2KM
            for idx in np.arange(1, dnew.shape[0]+1):
                tps[idx-1] = simpson(np.sqrt( np.power(vs[:idx], -2.) - p*p) -
                                     np.sqrt( np.power(vp[:idx], -2.) - p*p), dnew[:idx])
            #end for


            t.trim(starttime=t.stats.onset, endtime=t.stats.endtime)
            ipt.trim(starttime=ipt.stats.onset, endtime=ipt.stats.endtime)

            t2aio = interp1d(t.times(), t.data)
            t2irio = interp1d(ipt.times(), np.real(ipt.data))
            t2icio = interp1d(ipt.times(), np.imag(ipt.data))

            if(point_data_array is None):
                point_data_array = np.zeros((len(self._p_traces), NZ, 7))
            # end if

            point_data_array[it, :, 0] = np.squeeze(t2xio(d2tio(dnew)))
            point_data_array[it, :, 1] = np.squeeze(t2yio(d2tio(dnew)))
            point_data_array[it, :, 2] = np.squeeze(t2zio(d2tio(dnew)))
            point_data_array[it, :, 3] = np.squeeze(t2dio(d2tio(dnew)))

            point_data_array[it, :, 4] = np.squeeze(t2aio(tps))
            point_data_array[it, :, 5] = np.squeeze(t2irio(tps))
            point_data_array[it, :, 6] = np.squeeze(t2icio(tps))
        # end for
        #pbar.close()

        # serialize writing of rf_streams to disk
        for irank in np.arange(self._nproc):
            if(irank == self._rank):
                    if(len(point_data_array)):
                        point_data_array = point_data_array.reshape(-1, point_data_array.shape[-1])

                        hf = h5py.File(self._output_file, 'a')
                        dset = hf.create_dataset("%d"%(irank), point_data_array.shape, dtype=point_data_array.dtype)
                        dset[:] = point_data_array

                        for skey in self._stations.keys():
                            hf.attrs[skey] = self._stations[skey]
                        # end for

                        # add common metadata on rank 0
                        if(self._rank == 0):
                            hf.attrs['earth_radius'] = self._earth_radius
                        # end for

                        hf.close()
                    # end if
                # end if
            # end if
            self._comm.Barrier()
        # end for
    # end func
# end class


@click.command()
@click.option('--start-lat-lon', type=(float, float), required=True, help='Start latitude, longitude (degrees)')
@click.option('--azimuth', type=float, required=True, help='document me')
@click.option('--dimensions', type=(float, float, float), required=True,
              help='(length, width, depth) of the volume, in km')
@click.option('--num-nodes', type=(int, int, int), required=True,
              help='Number of nodes in each dimension of the volume')
@click.option('--debug', type=bool, is_flag=True, default=False)
@click.argument('rf-h5-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output_file', type=click.Path(exists=False), required=True)
def main(rf_h5_file, output_file, start_lat_lon, azimuth, dimensions, num_nodes, debug):
    """Perform 3D migration of RFs to volumetric space, stacking RF amplitudes in each cell.

    Example usage:
        python rf_3dmigrate.py --start-lat-lon -17.4 132.9 --azimuth 80 --dimensions 1000 450 75 \
            --num-cells 100 45 375 /g/data/ha3/am7399/shared/OA-ZRT-R-cleaned.h5 /g/data/ha3/am7399/shared/OA_piercing

    The script produces text data files which are converted to visualization using experimental
    ipython notebook `sandbox/plot_3dmigrate.ipynb`.

    :param rf_h5_file: Source file containing receiver functions
    :type rf_h5_file: str or Path
    :param output_folder: Folder in which to output results
    :type output_folder: str or Path
    """

    # Make sure geographiclib, needed for ray-tracing, exists
    if(not geodetics.HAS_GEOGRAPHICLIB):
        log.error("Python package 'geographiclib' not found. Aborting..")
        exit(0)
    # end if

    """
    g = Geometry(start_lat_lon=start_lat_lon, azimuth=azimuth,
                 lengthkm=dimensions[0], nx=num_nodes[0],
                 widthkm=dimensions[1], ny=num_nodes[1],
                 depthkm=dimensions[2], nz=num_nodes[2], debug=debug)
    """

    m = Migrate(rf_filename=rf_h5_file, debug=False, output_file=output_file)

    #m.execute()
# end

if __name__ == "__main__":
    # call main function
    main()  # pylint: disable=no-value-for-parameter
