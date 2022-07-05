"""
Description:
    Implements migration algorithm as described in Frassetto et al. (2010):
      Improved imaging with phase-weighted common conversion point stacks
      of receiver functions (GJI)

References:

CreationDate:   03/09/21
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     03/09/21   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from collections import defaultdict
from seismic.receiver_fn import rf_util, rf_corrections
import numpy as np
from rf import read_rf, RFStream
from rf.util import DEG2KM
from obspy.taup import TauPyModel
import h5py

from scipy.interpolate import splev, splrep, sproot, interp1d
from scipy.integrate import simps as simpson
from scipy.signal import hilbert
from mpi4py import MPI
from scipy.spatial import cKDTree
from pyproj import Geod
from shapely.geometry import Point, LineString, Polygon
import gdal
from osgeo.gdalconst import *
from affine import Affine
import struct
from seismic.stream_io import get_obspyh5_index

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

class Gravity:
    def __init__(self, gravity_grid_fn):
        self._fn = gravity_grid_fn
        self._ds = None
        try:
            self._ds = gdal.Open(self._fn, gdal.GA_ReadOnly)
            self._band = self._ds.GetRasterBand(1)
        except Exception as e:
            print(str(e))
            assert 0, 'Failed to load gravity grid. Aborting..'
        # end try

        self._gt = self._ds.GetGeoTransform()
        self._affine_forward_transform = Affine.from_gdal(*self._gt)
        self._affine_reverse_transform = ~(self._affine_forward_transform)
        # end func

    def _readPixel(self, px, py):
        def translateFormat(pt):
            fmttypes = {
                GDT_Byte: 'B',
                GDT_Int16: 'h',
                GDT_UInt16: 'H',
                GDT_Int32: 'i',
                GDT_UInt32: 'I',
                GDT_Float32: 'f',
                GDT_Float64: 'd'
            }
            return fmttypes.get(pt, 'x')

        # end func

        structval = self._band.ReadRaster(int(px), int(py), 1, 1, buf_type=self._band.DataType)
        val = struct.unpack(translateFormat(self._band.DataType), structval)[0]

        if val == self._band.GetNoDataValue():
            val = np.nan
        # end if

        return val
    # end func

    def query(self, lon, lat):
        # Convert geographic co-ordinates to pixel co-ordinates
        ipx, ipy = self._affine_reverse_transform * (lon, lat)
        ipx, ipy = int(ipx + 0.5), int(ipy + 0.5)

        return self._readPixel(ipx, ipy)
    # end func
# end class

class Migrate:
    def __init__(self, rf_filename, min_slope_ratio=-1, earth_radius=6371, logger=None):
        self._rf_filename = rf_filename
        self._min_slope_ratio = min_slope_ratio
        self._earth_radius = earth_radius #km
        self._logger = logger
        # Initialize MPI
        self._comm = MPI.COMM_WORLD
        self._nproc = self._comm.Get_size()
        self._rank = self._comm.Get_rank()

        self._p_traces = None # primary RF traces
        self._iphase_traces = None # instantaneous phase traces
        self._stations = defaultdict(list)
    # end func

    def process_streams(self, output_file, fmin=None, fmax=None, primary_component='R', model='iasp91', phase='P'):
        proc_hkeys = None
        if(self._rank == 0):
            hkeys = get_obspyh5_index(self._rf_filename, seeds_only=True)
            assert len(hkeys), 'No hdf5-groups found in file {}. Aborting..'.format(self._rf_filename)

            # split workload
            proc_hkeys = rf_util.split_list(hkeys, self._nproc)
        # end if
        # broadcast workload to all procs
        proc_hkeys = self._comm.bcast(proc_hkeys, root=0)

        # read local traces
        self._p_traces = RFStream()
        for hkey in proc_hkeys[self._rank]:
            if(self._logger): self._logger.info('rank {}: loading {}..'.format(self._rank, hkey))

            traces = read_rf(self._rf_filename, format='h5', group='waveforms/%s'%(hkey))

            # select primary component only
            p_traces = traces.select(component=primary_component)
            assert len(p_traces), 'No {} component found in RF stream for {}. Aborting..'.format(primary_component,
                                                                                                 hkey)
            if (self._min_slope_ratio > 0):
                before = len(p_traces)
                p_traces = RFStream([tr for tr in p_traces \
                                     if tr.stats.slope_ratio >= self._min_slope_ratio])
                after = len(p_traces)

                self._logger.info('rank {}: {} ({}/{}) traces '
                                  'dropped with min-slope-ratio filter..'.format(self._rank, hkey,
                                                                                 before-after, before))
            # end if

            if(len(p_traces) == 0):
                self._logger.warn('rank {}: {}: No traces left to process..'.format(self._rank, hkey))
                continue
            # end if

            has_reverberations = rf_corrections.has_reverberations(p_traces)
            if(has_reverberations):
                if (self._logger):
                    self._logger.info('rank {}: {}: removing reverberations..'.format(self._rank, hkey))
                # end if
                p_traces = rf_corrections.apply_reverberation_filter(p_traces)
            # end if

            if(fmin and fmax):
                self._logger.info('rank {}: {}: applying a bandpass filter..'.format(self._rank, hkey))
                p_traces.filter(type='bandpass', freqmin=fmin,
                                freqmax=fmax, corners=1, zerophase=True)
            elif(fmin):
                self._logger.info('rank {}: {}: applying a highpass filter..'.format(self._rank, hkey))
                p_traces.filter(type='highpass', freq=fmin,
                                corners=1, zerophase=True)
            elif(fmax):
                self._logger.info('rank {}: {}: applying a lowpass filter..'.format(self._rank, hkey))
                p_traces.filter(type='lowpass', freq=fmax,
                                corners=1, zerophase=True)
            # end if

            self._p_traces += p_traces
            self._stations[hkey] = [p_traces[0].stats.station_longitude,
                                    p_traces[0].stats.station_latitude,
                                    p_traces[0].stats.station_elevation,
                                    has_reverberations]
        # end for

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

        point_data_array = None
        for it, (t, ipt) in enumerate(zip(self._p_traces, self._iphase_traces)):
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

            #===============================================================
            # Create time-to-[depth,x,y,z] interpolants
            #===============================================================
            t2dio  = interp1d(times, depths)
            t2xio  = interp1d(times, xyzs[:,0])
            t2yio  = interp1d(times, xyzs[:,1])
            t2zio  = interp1d(times, xyzs[:,2])

            #===============================================================
            # Create a depth-to-time interpolant as expelined here:
            # https://stackoverflow.com/questions/47275957/invert-interpolation-to-give-the-variable-associated-with-a-desired-interpolatio
            #===============================================================
            tck = splrep(times, depths, k=3, s=0)

            NZ = 1500
            ZMAX = 150 # km
            dnew = np.linspace(1e-5, ZMAX, NZ)
            tnew = np.zeros(dnew.shape)
            for idx in np.arange(dnew.shape[0]):
                tck_mod = (tck[0], tck[1] - dnew[idx], tck[2])
                tnew[idx] = sproot(tck_mod)[0]
            #end for

            d2tio = interp1d(np.concatenate(([0], dnew)), np.concatenate(([0], tnew)))
            vp = vmodel.evaluate_above(dnew, 'p')
            vs = vmodel.evaluate_above(dnew, 's')

            #========================================================================
            # Create an interpolant for the delay time (t_Ps) of a non-vertically
            # incident Ps wave as follows:
            #          0
            # tps(d) = ∫ sqrt(Vs^-2 - p^2) - sqrt(Vp^-2 - p^2) dd
            #          d
            # Ref: https://ds.iris.edu/media/workshop/2013/01/advanced-studies-institute-on-seismological-research/files/Frassetto_RF-Stack_ASI.pdf
            # Note that the later phases (PpPs and PsPs+PpSs) are currently ignored.
            #========================================================================
            tps = np.zeros(dnew.shape)
            p = t.stats.slowness / DEG2KM
            for idx in np.arange(1, dnew.shape[0]+1):
                tps[idx-1] = simpson(np.sqrt( np.power(vs[:idx], -2.) - p*p) -
                                     np.sqrt( np.power(vp[:idx], -2.) - p*p), dnew[:idx])
            #end for


            t.trim(starttime=t.stats.onset, endtime=t.stats.endtime)
            ipt.trim(starttime=ipt.stats.onset, endtime=ipt.stats.endtime)

            #=========================================================================
            # Create time-to-data interpolants for RF amplitude and real and imaginary
            # components of the analytic representation
            #=========================================================================
            t2aio = interp1d(t.times(), t.data, bounds_error=False, fill_value=0)
            t2irio = interp1d(ipt.times(), np.real(ipt.data), bounds_error=False, fill_value=0)
            t2icio = interp1d(ipt.times(), np.imag(ipt.data), bounds_error=False, fill_value=0)

            # array to store values for all station traces to be written to an HDF5 file
            if(point_data_array is None):
                point_data_array = np.zeros((len(self._p_traces), NZ, 7))
            # end if

            #=========================================================================
            # Populate array with values for x, y, z and depth, going from:
            # depth -> time -> x, etc.
            #=========================================================================
            point_data_array[it, :, 0] = np.squeeze(t2xio(d2tio(dnew)))
            point_data_array[it, :, 1] = np.squeeze(t2yio(d2tio(dnew)))
            point_data_array[it, :, 2] = np.squeeze(t2zio(d2tio(dnew)))
            point_data_array[it, :, 3] = np.squeeze(t2dio(d2tio(dnew)))

            #==========================================================================
            # Populate array with RF trace amplitudes and real and imaginary components
            # of the analytic representation at t_Ps delay times
            #==========================================================================
            point_data_array[it, :, 4] = np.squeeze(t2aio(tps))
            point_data_array[it, :, 5] = np.squeeze(t2irio(tps))
            point_data_array[it, :, 6] = np.squeeze(t2icio(tps))
        # end for

        # serialize writing of rf_streams to disk
        for irank in np.arange(self._nproc):
            if(irank == self._rank):
                    if(len(point_data_array)):
                        point_data_array = point_data_array.reshape(-1, point_data_array.shape[-1])

                        if(self._logger): self._logger.info('rank {}: saving migrated volume.. '.format(self._rank))

                        hf = h5py.File(output_file, 'a')
                        dset = hf.create_dataset("%d"%(irank), point_data_array.shape, dtype=point_data_array.dtype)
                        dset[:] = point_data_array

                        for skey in self._stations.keys():
                            dset.attrs[skey] = self._stations[skey]
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
        if(self._rank == 0):
            if (self._logger): self._logger.info('Done..')
        # end if
    # end func
# end class

class CCPVolume():
    def __init__(self, fn):
        self._fn = fn
        self._meta = defaultdict(list)
        self._earth_radius = None
        self._data = None
        self._tree = None

        hf = None
        try:
            hf = h5py.File(self._fn, 'r')
        except Exception as e:
            print(str(e))

            assert 0, 'Failed to load {}. Aborting..'.format(self._fn)
        # end try

        # read metadata
        for k in hf.attrs.keys():
            if (k == 'earth_radius'):
                self._earth_radius = hf.attrs[k]
        # end for

        # read station metadata
        for dk in hf.keys():
            for sk in hf[dk].attrs.keys():
                self._meta[sk] = hf[dk].attrs[sk] # station -> lon, lat, elevation, hasReverberations
            # end for
        # end for

        # read ccp volume
        self._data = []
        for k in hf.keys():
            self._data.append(np.array(hf[k]))
        # end for
        self._data = np.vstack(self._data)

        self._tree = cKDTree(self._data[:, :3], balanced_tree=False)

        hf.close()
    # end func
# end class

class CCP_VerticalProfile():
    def __init__(self, ccpVolume, profile_start, profile_end, dx=5, dz=0.5, max_depth=150,
                 swath_width=40, ds=10, extend=50, cell_radius=20, idw_exponent=2,
                 pw_exponent=1, max_station_dist=10):
        self._ccpVolume = ccpVolume

        self._profile_start = profile_start
        self._profile_end = profile_end
        self._lon1 = None
        self._lat1 = None
        self._lon2 = None
        self._lat2 = None
        self._r1 = None
        self._r2 = None
        self._dx = dx
        self._dz = dz
        self._max_depth = max_depth
        self._swath_width = swath_width
        self._ds = ds
        self._extend = extend
        self._cell_radius = cell_radius
        self._idw_exponent = idw_exponent
        self._pw_exponent = pw_exponent
        self._max_station_dist = max_station_dist
        self._profileLength = None
        self._nLateralNodes = None
        self._nDepthNodes = None
        self._nSwathNodes = None
        self._lateralNodes = None
        self._depthNodes = None

        self._grid = None
        self._gx = None
        self._gd = None
        self._gs = None
        self._gmeta = None
        self._grid_vals = None

        # fetch coordinates from CCP-volume
        def get_location(loc):
            lon = lat = None
            try:
                if(type(loc) == str):
                    try:
                        lon, lat = self._ccpVolume._meta[loc][:2]
                    except Exception as e:
                        print(str(e))
                        assert 0, 'Station not found. Aborting..'
                    # end try
                else:
                    lon, lat = loc[0], loc[1]
                # end if
            except Exception as e:
                print(str(e))
                assert 0, 'Error encountered processing profile start/end location. Aborting..'
            # end try

            return lon, lat
        # end func

        lon1, lat1 = get_location(self._profile_start)
        lon2, lat2 = get_location(self._profile_end)

        self._generateGrid(lon1, lat1, lon2, lat2)
        self._generateGridMeta()
        self._interpolate()
    # end func

    def _generateGridMeta(self):
        # generate profile-meta
        sta_list = []
        coords = []
        meta = self._ccpVolume._meta
        ER = self._ccpVolume._earth_radius

        # get all station coordinates in CCP volume
        for k in meta.keys():
            if (len(meta[k])):
                sta_list.append(k)
                coords.append(meta[k][:2])
            # end if
        # end for

        # create a kd-Tree of station coordinates
        coords = np.array(coords)
        coords_xyz = rtp2xyz(np.ones(coords.shape[0]) * ER,
                             np.radians(90 - coords[:, 1]),
                             np.radians(coords[:, 0]))
        ct = cKDTree(coords_xyz)

        # for all nodes in the profile, get distance and index for
        # closest station
        nodes_xyz = rtp2xyz(np.ones(self._lateralNodes.shape[0]) * ER,
                            np.radians(90 - self._lateralNodes[:, 1]),
                            np.radians(self._lateralNodes[:, 0]))
        dlist, llist = ct.query(nodes_xyz)

        profile_meta = defaultdict(list)
        for i, (d, l) in enumerate(zip(dlist, llist)):
            profile_meta[sta_list[l]].append([d, l, i])
        # end for

        # For each station, filter out all but the closest grid-node
        self._g_meta = {}
        for k in profile_meta.keys():
            arr = np.array(profile_meta[k])
            # print([k,arr])
            if (np.min(arr[:, 0]) > self._max_station_dist): continue
            sta_idx = np.int_(arr[np.argmin(arr[:, 0]), 1])
            node_idx = np.int_(arr[np.argmin(arr[:, 0]), 2])
            # print([gx[node_idx], sta_list[sta_idx]])
            self._g_meta[sta_list[sta_idx]] = {'distance': self._gx[node_idx], 'corrected': meta[k][-1]}
        # end for
    # end func

    def _generateGrid(self, lon1, lat1, lon2, lat2):
        # Compute lateral nodes
        geod = Geod(ellps="WGS84")

        # extend profile
        az, baz, _ = geod.inv(lon1, lat1, lon2, lat2)
        self._lon1, self._lat1, _ = geod.fwd(lon1, lat1, baz, self._extend * 1e3)
        self._lon2, self._lat2, _ = geod.fwd(lon2, lat2, az, self._extend * 1e3)

        self._r1 = self._ccpVolume._earth_radius - self._max_depth
        self._r2 = self._ccpVolume._earth_radius

        self._profileLength = geod.geometry_length(LineString([Point(self._lon1, self._lat1),
                                                               Point(self._lon2, self._lat2)])) / 1e3
        self._nLateralNodes = np.int_(np.ceil(self._profileLength / self._dx)) + 1
        self._lateralNodes = np.vstack([np.array([self._lon1, self._lat1]),
                                        np.array(geod.npts(self._lon1, self._lat1,
                                                 self._lon2, self._lat2,
                                                 self._nLateralNodes-2)),
                                        np.array([self._lon2, self._lat2])])

        # swath nodes
        self._nSwathNodes = int((self._swath_width / 2. / self._ds) * 2 + 1)
        self._gs = np.linspace(-self._swath_width / 2., self._swath_width / 2., self._nSwathNodes)

        # depth nodes
        self._nDepthNodes = np.int_(np.ceil((self._r2 - self._r1) / self._dz)) + 1
        self._depthNodes = np.linspace(self._r1, self._r2, self._nDepthNodes)

        # Assemble lon-lat nodes for profile
        self._grid = np.zeros((self._nLateralNodes * self._nDepthNodes, 3))
        for i in np.arange(self._nDepthNodes):
            for j in np.arange(self._nLateralNodes):
                idx = i * self._nLateralNodes + j
                self._grid[idx, :] = np.array([self._depthNodes[i],
                                               self._lateralNodes[j, 0],
                                               self._lateralNodes[j, 1]])
            # end for
        # end for

        self._gx = np.linspace(0, self._profileLength, self._nLateralNodes)
        self._gd = np.linspace(self._r2 - self._r1, 0, self._nDepthNodes)

    # end func

    def _interpolate(self):
        ER = self._ccpVolume._earth_radius
        gxyz = rtp2xyz(self._grid[:, 0],
                       np.radians(90 - self._grid[:, 2]),
                       np.radians(self._grid[:, 1]))

        # iterate over swath node
        swath_grid_vals = []
        swath_grid_iphase_weights = []
        for s in self._gs:
            # compute swath-grid (profile nodes are simply shifted along swath)
            a = gxyz[0]
            b = gxyz[(self._gd.shape[0] - 1) * self._gx.shape[0]]
            c = gxyz[-1]
            ab = b - a
            ac = c - a
            normal = np.cross(ab, ac)
            normal /= np.linalg.norm(normal)

            sxyz = gxyz + normal * s

            # temporary arrays
            grid_vals = np.zeros(sxyz.shape[0])
            grid_iphase_weights = np.zeros(sxyz.shape[0], dtype=np.complex)

            p = self._idw_exponent
            r = self._cell_radius
            tree = self._ccpVolume._tree
            data = self._ccpVolume._data

            # iterate over nodes in each swath
            for i in np.arange(sxyz.shape[0]):

                indices = np.array(tree.query_ball_point(sxyz[i, :], r=r))

                if (len(indices) == 0): continue

                # filter out all but nodes within a disc of dz thickness
                indices = np.array(indices)
                indices = indices[np.fabs(data[indices, 3] -
                                          np.fabs((self._grid[i, 0] - ER))) < self._dz]

                if (len(indices) == 0): continue
                d = np.zeros(len(indices))

                # compute distance of kdtree nodes from current node in swath
                d[:] = np.sqrt(np.sum(np.power(sxyz[i] - data[indices, :3], 2), axis=1))

                # filter out nodes outside a cone, defined as current_radius = current depth;
                # this is done to avoid lateral smearing at shallow depths, where piercing
                # points are sparse, except immediately under stations
                indices = indices[d < np.fabs((self._grid[i, 0] - ER))]
                d = d[d < np.fabs((self._grid[i, 0] - ER))]

                if (len(indices) == 0): continue

                # compute IDW weights
                idwIndices = indices
                idw = np.zeros(d.shape)
                idw = 1. / np.power(d, p)

                # compute mean instantaneous phase weight
                grid_iphase_weights[i] = np.mean(data[idwIndices, 5] + 1j * data[idwIndices, 6])

                # compute grid values
                grid_vals[i] = np.sum(idw * data[idwIndices, 4]) / np.sum(idw)
            # end for

            swath_grid_iphase_weights.append(grid_iphase_weights)
            swath_grid_vals.append(grid_vals)
        # end for

        # sum grid values across swath
        v = np.reshape(np.array(swath_grid_vals),
                       (self._gs.shape[0], self._gd.shape[0], self._gx.shape[0]))
        pw = np.reshape(np.array(swath_grid_iphase_weights),
                       (self._gs.shape[0], self._gd.shape[0], self._gx.shape[0]))

        self._grid_vals = np.mean(v * np.power(np.abs(pw), self._pw_exponent), axis=0)
        #print(self._grid_vals)
    # end func

    def plot(self, ax, amp_min=-0.2, amp_max=0.2, gax=None, gravity=None):
        if(gax and gravity):
            gvs = np.zeros(self._gx.shape)
            for i in np.arange(self._nLateralNodes):
                gvs[i] = gravity.query(self._lateralNodes[i, 0], self._lateralNodes[i, 1])
            # end for

            gax.plot(self._gx, gvs)
            for pos in ['right', 'top']:
                gax.spines[pos].set_visible(False)
            # end for
        # end if

        tickstep_x = 50
        tickstep_y = 5

        ax.set_xlabel('Distance [km]', fontsize=12)
        ax.set_ylabel('Depth [km]', fontsize=12)
        ax.tick_params(direction="in", labelleft=True, labelright=True)

        ax.set_xticks(np.arange(0, np.max(self._gx) * 1.0001, tickstep_x))
        ax.set_yticks(np.arange(0, np.max(self._gd) * 1.0001, tickstep_y))

        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')

        ax.grid(color='#808080', which='major', linestyle=':', alpha=0.5)

        cs = ax.contourf(self._gx, self._gd,
                         self._grid_vals,
                         cmap='seismic',
                         levels=np.linspace(amp_min, amp_max, 100),
                         extend='both')
        ax.invert_yaxis()
        for k in self._g_meta.keys():
            px = self._g_meta[k]['distance']
            py = -4.
            ax.text(px, py, "{}{}".format('*' if self._g_meta[k]['corrected'] else '', k),
                    horizontalalignment='center',
                    verticalalignment='top', fontsize=9, backgroundcolor='#ffffffa0')
        # end for

        titleAx = None
        if(gax): titleAx = gax
        else: titleAx = ax

        titleAx.set_title('CCP Stack: {} --- {}'.format(str(self._profile_start), str(self._profile_end)),
                          fontdict={'fontsize': 14}, pad=40)
    # end func
# end class

class CCP_DepthProfile():
    def __init__(self, ccpVolume, net_sta_loc1, net_sta_loc2, depth, dx=5, dy=5, dz=1,
                 extend=50, cell_radius=40, idw_exponent=2, pw_exponent=1):
        """
        :param ccpVolume:
        :param profile_start: Either a NET.STA.LOC string or a tuple (lon, lat)
        :param profile_end: Either a NET.STA.LOC string or a tuple (lon, lat)
        :param depth:
        :param dx:
        :param dy:
        :param dz:
        :param extend:
        :param cell_radius:
        :param idw_exponent:
        :param pw_exponent:
        """
        self._ccpVolume = ccpVolume

        self._net_sta_loc1 = net_sta_loc1
        self._net_sta_loc2 = net_sta_loc2
        self._lon1 = None
        self._lat1 = None
        self._lon2 = None
        self._lat2 = None
        self._depth = depth
        self._dx = dx
        self._dy = dy
        self._dz = dz
        self._extend = extend
        self._cell_radius = cell_radius
        self._idw_exponent = idw_exponent
        self._pw_exponent = pw_exponent
        self._length_x = None
        self._length_y = None
        self._nx = None
        self._ny = None

        self._ynodes = None
        self._grid = None
        self._gx = None
        self._gy = None
        self._gs = None
        self._gmeta = None
        self._grid_vals = None

        # fetch coordinates from CCP-volume
        lon1 = lat1 = lon2 = lat2 = None
        try:
            lon1, lat1 = self._ccpVolume._meta[net_sta_loc1][:2]
            lon2, lat2 = self._ccpVolume._meta[net_sta_loc2][:2]
        except Exception as e:
            print(str(e))
            assert 0, 'Station not found. Aborting..'
        # end try

        self._generateGrid(lon1, lat1, lon2, lat2)
        self._interpolate()

    # end func

    def _generateGrid(self, lon1, lat1, lon2, lat2):
        geod = Geod(ellps="WGS84")

        # extend profile
        az, baz, _ = geod.inv(lon1, lat1, lon2, lat2)

        self._lon1, self._lat1, _ = geod.fwd(lon1, lat1, baz, self._extend * 1e3)
        self._lon2, self._lat2, _ = geod.fwd(lon2, lat2, az, self._extend * 1e3)

        self._length_x = geod.geometry_length(LineString([Point(self._lon1, self._lat1),
                                                          Point(self._lon2, self._lat1)])) / 1e3
        self._length_y = geod.geometry_length(LineString([Point(self._lon1, self._lat1),
                                                          Point(self._lon1, self._lat2)])) / 1e3

        self._nx = np.int_(np.ceil(self._length_x / self._dx)) + 1
        self._ny = np.int_(np.ceil(self._length_y / self._dy)) + 1

        self._ynodes = np.vstack([np.array([self._lon1, self._lat1]),
                                  np.array(geod.npts(self._lon1, self._lat1,
                                                     self._lon1, self._lat2,
                                                     self._ny - 2)),
                                  np.array([self._lon1, self._lat2])])

        r = self._ccpVolume._earth_radius - self._depth
        az, baz, _ = geod.inv(self._lon1, self._lat1, self._lon2, self._lat1)

        # Assemble 3D nodes
        self._grid = np.zeros((self._nx * self._ny, 3))
        for i in np.arange(self._ny):
            startLon, startLat = self._ynodes[i]
            endLon, endLat, _ = geod.fwd(startLon, startLat, az, self._length_x * 1e3)
            xnodes = np.vstack([np.array([startLon, startLat]),
                                np.array(geod.npts(startLon, startLat, endLon, endLat, self._nx - 2)),
                                np.array([endLon, endLat])])

            for j in np.arange(self._nx):
                idx = i * self._nx + j
                self._grid[idx, :] = np.array([r, xnodes[j, 0], xnodes[j, 1]])
            # end for
        # end for

        self._gx = np.linspace(0, self._length_x, self._nx)
        self._gy = np.linspace(0, self._length_y, self._ny)

        # Gather metadata
        p = Polygon([(self._lon1, self._lat1), (self._lon2, self._lat1),
                     (self._lon2, self._lat2), (self._lon1, self._lat2)])
        self._gmeta = {}
        for k in self._ccpVolume._meta.keys():
            attrs = self._ccpVolume._meta[k]

            site = Point(attrs[0], attrs[1])
            if (p.contains(site)):
                distx = geod.geometry_length(LineString([Point(self._lon1, self._lat1),
                                                         Point(attrs[0], self._lat1)])) / 1e3
                disty = geod.geometry_length(LineString([Point(self._lon1, self._lat1),
                                                         Point(self._lon1, attrs[1])])) / 1e3

                self._gmeta[k] = {'distx': distx, 'disty': disty, 'corrected': attrs[-1]}
            # end if
        # end for

    # end func

    def _interpolate(self):
        ER = self._ccpVolume._earth_radius
        gxyz = rtp2xyz(self._grid[:, 0],
                       np.radians(90 - self._grid[:, 2]),
                       np.radians(self._grid[:, 1]))

        grid_vals = np.zeros(gxyz.shape[0])
        grid_iphase_weights = np.zeros(len(gxyz), dtype=np.complex)

        p = self._idw_exponent
        r = self._cell_radius
        tree = self._ccpVolume._tree
        data = self._ccpVolume._data

        for i in np.arange(gxyz.shape[0]):
            indices = np.array(tree.query_ball_point(gxyz[i, :], r=r))

            if (len(indices) == 0): continue

            indices = np.array(indices)

            indices = indices[np.fabs(data[indices, 3] - (ER - self._grid[i, 0])) < self._dz]
            if (len(indices) == 0): continue

            d = np.zeros(len(indices))
            d[:] = np.sqrt(np.sum(np.power(gxyz[i] - data[indices, :3], 2), axis=1))

            idwIndices = indices
            idw = np.zeros(d.shape)
            idw = 1. / np.power(d, p)

            grid_iphase_weights[i] = np.mean(data[idwIndices, 5] + 1j * data[idwIndices, 6])
            grid_vals[i] = np.sum(idw * data[idwIndices, 4]) / np.sum(idw) * \
                           np.power(np.abs(grid_iphase_weights[i]), self._pw_exponent)
        # end for

        self._grid_vals = grid_vals

    # end func

    def plot(self, ax, amp_min=-0.2, amp_max=0.2):
        ax.set_xlabel('Distance [km]', fontsize=12)
        ax.set_ylabel('Distance [km]', fontsize=12)

        cs = ax.contourf(self._gx, self._gy,
                         np.reshape(self._grid_vals, (self._gy.shape[0], self._gx.shape[0])),
                         cmap='seismic', levels=np.linspace(amp_min, amp_max, 100))
        for c in cs.collections:
            c.set_rasterized(True)
        # end for

        for k in self._gmeta.keys():
            distx = self._gmeta[k]['distx']
            disty = self._gmeta[k]['disty']
            ax.text(distx, disty, "{}{}".format('*' if self._gmeta[k]['corrected'] else '', k),
                    horizontalalignment='right', alpha=0.2,
                    bbox=dict(facecolor='none', alpha=0.4, edgecolor='none'),
                    fontsize=4, backgroundcolor='#ffffffa0')
            ax.scatter(distx, disty, s=0.1, marker='x', c='k')
        # end for

        ax.set_title('CCP Depth profile: {} km'.format(self._depth),
                     fontdict={'fontsize': 14}, pad=40)

        #print(np.reshape(self._grid_vals, (self._gy.shape[0], self._gx.shape[0])))
    # end func
# end class
