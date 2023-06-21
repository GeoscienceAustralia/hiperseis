"""
Description:
    A collection of segy-related utilities and functions

References:

CreationDate:   19/06/23
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     19/06/23   RH
"""

from scipy.interpolate import interp1d
import numpy as np
import pyproj

from obspy.io.segy import segy

class DepthMigratedSegy:
    def __init__(self, segy_fn,
                 depth_zero_km,
                 max_depth_km,
                 coords_file=None,
                 epsg_code=28353,
                 pick_every=1):
        '''
        Class for reading 2D segy files. This class assumes that the traces are
        CDP-ordered.

        :param segy_fn: segy file name
        :param depth_zero_m: depth of first sample in kilometres
        :param max_depth_m: maximum depth in kilometres
        :param coords_file: comma-separated text file containing x and y coordiates to be read
                            in instead of using coordinates available in segy file.
        :param epsg_code: epsg code for projection used (default is GDA94 / MGA zone 53)
        :param pick_every: by default, every trace is read from the segy file, but
                           pick_every>1 allows skipping traces.
        '''
        self.sfn = segy_fn
        self.coords_file = coords_file
        self.depth_zero_km = depth_zero_km
        self.max_depth_km = max_depth_km
        self.epsg_code = epsg_code
        self.sf = segy._read_segy(self.sfn)
        self.pick_every = pick_every
        # transformer to go from segy projection to wgs84
        self.transformer = pyproj.Transformer.from_crs(self.epsg_code, 4326)

        # Read traces
        self.samples = []  # trace samples
        self.xs = []  # x coordinates
        self.ys = []  # y coordinates
        self.ns = []  # num samples
        self.cdps = []  # ensemble number (cdps)
        self.si = []  # sample interval
        self.ntraces = 0
        for itr, tr in enumerate(self.sf.traces):
            if (itr % self.pick_every == 0):
                self.samples.append(tr.data)
                self.xs.append(tr.header.source_coordinate_x)
                self.ys.append(tr.header.source_coordinate_y)
                self.ns.append(tr.header.number_of_samples_in_this_trace)
                self.si.append(tr.header.sample_interval_in_ms_for_this_trace / 1e6)  # convert from micro seconds to s
                self.cdps.append(tr.header.ensemble_number)
                self.ntraces += 1
        # end for

        if(self.coords_file):
            coords = np.loadtxt(self.coords_file, delimiter=',')
            self.xs[:] = coords[:, 0]
            self.ys[:] = coords[:, 1]
        # end if

        self.samples = np.array(self.samples).T  # shape (nz, nt)
        self.ns = np.array(self.ns)
        assert np.min(self.ns) == np.max(self.ns), \
            'Sample-count mismatch found. Aborting..'
        self.si = np.array(self.si)
        self.xs = np.array(self.xs)
        self.ys = np.array(self.ys)
        self.cdps = np.array(self.cdps)
        self.station_spacing = np.sqrt((self.xs[:-1] - self.xs[1:]) ** 2 +
                                       (self.ys[:-1] - self.ys[1:]) ** 2)
        self.station_spacing = np.concatenate([[0], self.station_spacing])

        # use station-spacing to weed out traces with incorrect coordinates
        median_spacing = np.median(self.station_spacing)
        bad_indices = np.fabs(self.station_spacing - median_spacing) > 0.1*median_spacing
        print('Warning: removing {} traces (~{:.3f}%) with '
              'bad coordinates..'.format(np.sum(bad_indices),
                                         np.sum(bad_indices)/self.ntraces*100))

        self.ns = self.ns[~bad_indices]
        self.si = self.si[~bad_indices]
        self.xs = self.xs[~bad_indices]
        self.ys = self.ys[~bad_indices]
        self.cdps = self.cdps[~bad_indices]
        self.station_spacing= self.station_spacing[~bad_indices]
        self.ntraces -= np.sum(bad_indices)
        self.samples = self.samples[:, ~bad_indices]

        # compute euclidean distance along profile
        self.ds = np.array([np.sum(self.station_spacing[:i])
                            for i in range(len(self.station_spacing))]) / 1e3  # km
        self.ds -= np.min(self.ds)  # distance along profile starts from 0

        self.times = np.linspace(0, (np.max(self.ns) - 1) * np.max(self.si),
                                 np.max(self.ns))
        # depths
        self.zs = np.linspace(self.depth_zero_km, self.max_depth_km,
                              np.max(self.ns))

        # convert x and y coordinates to lons and lats
        self.lats, self.lons = self.transformer.transform(self.xs, self.ys)
    # end func

    def getAttribute(self, key, d):
        '''
        Returns attribute value -- for the given key -- of the closest trace at a given
        distance along the seismic profile.

        :param key: attribute key; see key-value pairs below:
                    'samples' : trace samples
                    'x': x-coordinate of trace
                    'y': y-coordinate of trace
                    'lon': longitude of trace
                    'lat': latitude of trace
                    'cdp': CDP of trace
        :param d: distance along profile in km
        :return: attribute value of trace samples for a trace at a given distance along the
                 seismic profile.
        '''

        if (key not in ['samples', 'x', 'y', 'lon', 'lat', 'cdp']):
            assert 0, "Invalid key; should be one of \
                     ['samples', 'x', 'y', 'lon', 'lat', 'cdp']"
        if (d <= np.max(self.ds)):
            idx = np.argmin(np.fabs(self.ds - d))

            if (key == 'samples'):
                return self.samples[:, idx]
            elif (key == 'x'):
                return self.xs[idx]
            elif (key == 'y'):
                return self.ys[idx]
            elif (key == 'lon'):
                return self.lons[idx]
            elif (key == 'lat'):
                return self.lats[idx]
            elif (key == 'cdp'):
                return self._cdps[idx]
            else:
                raise ValueError('Attribute "%s" not found' % key)
            # end if
        else:
            return None
        # end if

    # end func

    def getProfile(self, ndepths=-1, nsteps=-1):
        '''
        Rreturns the seismic profile along with trace coordinates.

        :param ndepths: number of depth values to return for each trace,
                        starting from 0 to 'max_depth_km'. When set to -1, ndepth
                        is internally reset to the number of samples in a trace
        :param nsteps: number of steps along the profile. When set to -1, all traces
                       along the profile are returned, otherwise interpolated traces
                       at nsteps locations along the profile are returned
        :return: lonlat_list: 2D array of trace coordinates along the seismic line
                              of shape (ntraces, 2)
                 depth_list: depths for each sample of shape (ndepths)
                 distances: distances along the profile
                 samples: A 2D array of amplitudes of shape (ndepths, ntraces)
        '''

        if (ndepths < 0): ndepths = np.max(self.ns)
        if (nsteps < 0): nsteps = self.ntraces

        zs = np.linspace(self.depth_zero_km, self.max_depth_km, ndepths)
        ds = np.linspace(0, self.ds[-1], nsteps)

        result_samples = np.zeros((ndepths, nsteps))
        lons = np.zeros(nsteps)
        lats = np.zeros(nsteps)
        for i, d in enumerate(ds):
            samples = self.getAttribute('samples', d)
            lon = self.getAttribute('lon', d)
            lat = self.getAttribute('lat', d)

            io = interp1d(self.zs, samples,
                          bounds_error=False, fill_value=0)
            result_samples[:, i] = io(zs)
            lons[i] = lon
            lats[i] = lat
        # end for

        return lons, lats, zs, ds, result_samples
    # end func
# end class
