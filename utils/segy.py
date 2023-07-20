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
                 pick_every=1,
                 filter_coordinates=True,
                 filter_coordinates_factor=0.1,
                 flip=False):
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
        :param filter_coordinates: drop bad coordinates where station spacing exceeds
               filter_coordinates_factor * median_station_spacing
        :param filter_coordinates_factor: factor applied to median station spacing, while
               identifying problematic station-coordinates
        :param flip: Flip horizontally
        '''
        self.sfn = segy_fn
        self.coords_file = coords_file
        self.depth_zero_km = depth_zero_km
        self.max_depth_km = max_depth_km
        self.epsg_code = epsg_code
        self.sf = segy._read_segy(self.sfn)
        self.pick_every = pick_every
        self.filter_coordinates = filter_coordinates
        self.filter_coordinates_factor = filter_coordinates_factor
        self.flip = flip
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
            assert len(coords) == len(self.xs), 'Unexpected number of coordinates in {}. ' \
                                                'Aborting..'.format(self.coords_file)
            self.xs = coords[:, 0]
            self.ys = coords[:, 1]
        # end if

        self.samples = np.array(self.samples).T  # shape (nz, nt)

        # flip horizontally if required
        if(self.flip):
            self.xs = self.xs[::-1]
            self.ys = self.ys[::-1]
            self.samples = self.samples[:, ::-1]
        # end if

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
        bad_indices = np.zeros(self.station_spacing.shape, dtype='?')
        if(self.filter_coordinates):
            median_spacing = np.median(self.station_spacing)
            bad_indices = np.fabs(self.station_spacing - median_spacing) > \
                          self.filter_coordinates_factor * median_spacing
            print('Warning: removing {} traces (~{:.3f}%) with '
                  'bad coordinates..'.format(np.sum(bad_indices),
                                         np.sum(bad_indices)/self.ntraces*100))
        # end if

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

    def augment(self, other, prepend=False):
        """
        @param other: another segy line to append/prepend to this object
        @param prepend: prepend line, instead of appending it by default
        """
        # check compatibility
        if ((np.median(self.si) != np.median(other.si)) or
                (self.depth_zero_km != other.depth_zero_km) or
                (self.max_depth_km != other.max_depth_km)):
            print('DepthMigratedSegy intances are not compatible. '
                  'Aborting prepend operation..')
            return
        else:
            this_max, other_max = np.median(np.fabs(self.samples)), \
                                    np.median(np.fabs(other.samples))
            scale = this_max / other_max

            # append/prepend coordinates and data in place
            if (prepend):
                self.xs = np.hstack([other.xs, self.xs])
                self.ys = np.hstack([other.ys, self.ys])
                self.station_spacing = np.hstack([other.station_spacing, self.station_spacing])

                result = np.zeros((self.samples.shape[0],
                                   self.samples.shape[1] + other.samples.shape[1]))

                result[:, other.samples.shape[1]:] = self.samples[:, :]
                for i in np.arange(other.samples.shape[0]):
                    if (i >= self.samples.shape[0]): break
                    result[i, :other.samples.shape[1]] = other.samples[i, :] * scale
                # end for
                self.samples = result
            else:
                self.xs = np.hstack([self.xs, other.xs])
                self.ys = np.hstack([self.ys, other.ys])
                self.station_spacing = np.hstack([self.station_spacing, other.station_spacing])

                result = np.zeros((self.samples.shape[0],
                                   self.samples.shape[1] + other.samples.shape[1]))

                result[:, :self.samples.shape[1]] = self.samples[:, :]
                for i in np.arange(other.samples.shape[0]):
                    if (i >= self.samples.shape[0]): break
                    result[i, self.samples.shape[1]:] = other.samples[i, :] * scale
                # end for
                self.samples = result
            # end if

            # recompute euclidean distance along profile
            self.ds = np.array([np.sum(self.station_spacing[:i])
                                for i in range(len(self.station_spacing))]) / 1e3  # km

            # convert x and y coordinates to lons and lats
            self.lats, self.lons = self.transformer.transform(self.xs, self.ys)
        # end if
    # end func
# end class
