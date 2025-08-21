import numpy as np
import pandas as pd
from pandas.core.series import Series
from scipy.spatial import cKDTree
from pyproj import Geod
from obspy.core import UTCDateTime
from collections import defaultdict
from obspy.geodetics.base import degrees2kilometers
from seismic.misc import rtp2xyz
from obspy.core.event import Event, Origin, Magnitude, Catalog
from obspy.core.utcdatetime import UTCDateTime
from itertools import product
from typing import Tuple
from tqdm import tqdm

def azimuth_difference(a, b, p, ellipse='WGS84'):
    geod = Geod(ellps=ellipse)

    # ab    Geodesic from A to B.
    # az_ab Azimuth of geodesic from A to B.
    az_ab, _, _ = geod.inv(a[0], a[1], b[0], b[1])

    # ap    Geodesic from A to P.
    # az_ap Azimuth of geodesic from A to P.
    az_ap, _, d_ap = geod.inv(a[0], a[1], p[0], p[1])

    # bp    Geodesic from B to P.
    # az_bp Azimuth of geodesic from B to P.
    az_bp, _, d_bp = geod.inv(b[0], b[1], p[0], p[1])

    if(d_ap > d_bp):
        return np.min(np.array([np.fabs(az_ab - az_ap), np.fabs(az_ab-az_bp)]))
    else:
        return azimuth_difference(b, a, p)
    # end if
# end func

class GCMTCatalog:
    def __init__(self, source, ellipse='WGS84'):
        def read_gcmt_catalog(fn):
            """
            @param fn: GCMT catalog file name in text format. The expected columns (space-separated) are:
            lonp lon lat dep mrr mtt mff mrt mrf mtf exp EventOrigintim DC CLVD VOL Mw str1 dip1 rak1 \
            str2 dip2 rak2 plunP azP plunT azT Hlonp Hlon Hlat Hdep Ctime Chdur MwG base Bst Bco Bmp Bper \
            Sst Sco Smp Sper Mst Mco Mmp Mper
            @return: catalog in a pandas dataframe
            """
            def convert_to_timestamp(eotime):
                a = str(eotime)
                b = a[0:4] + '-' + a[4:6] + '-' + a[6:8] + 'T' + a[8:10] + ':' + a[10:12] + ':' + a[12:]
                return UTCDateTime(b).timestamp
            # end func

            cat = pd.read_csv(fn, header=[1], delimiter='\s+')
            cat['EventOrigintim'] = cat['EventOrigintim'].map(convert_to_timestamp)

            return cat
        # end func

        if (type(source) == pd.DataFrame):
            self.cat = source.copy()
        else:
            self.cat = read_gcmt_catalog(source)
        # end if
        self._initialize(ellipse=ellipse)
    # end func

    def _initialize(self, ellipse='WGS84'):
        self.EARTH_RADIUS_KM = 6371.

        self.geod = Geod(ellps=ellipse)
        self.tree = None

        # create kDTree for spatial queries
        r = np.ones(len(self.cat)) * self.EARTH_RADIUS_KM
        t = np.radians(90 - self.cat['lat'])
        p = np.radians(self.cat['lon'])

        xyz = rtp2xyz(r, t, p)
        self.tree = cKDTree(xyz)
    # end func

    def prune(self,
              time_range: Tuple[UTCDateTime, UTCDateTime] = None,
              lon=None, lat=None, distance_range: Tuple[float, float] = None,
              mag_range: Tuple[float, float] = None,
              depth_range: Tuple[float, float] = None,
              min_areal_separation_km=0):
        """
        @param time_range: start- and end-times to clip catalogue to
        @param lon: longitude to be used for restricting events to distance_range
        @param lat: latitude to be used for restricting events to distance_range
        @param distance_range: distance range to clip events to in degrees
        @param mag_range: magnitude Mw range to clip events to
        @param depth_range: depth range in km to clip events to
        @param min_areal_separation_km: areal extent over which proximal events that are within
                                        15 minutes of each other and for which magnitudes do not
                                        differ by more than 0.3 are dropped as being duplicates
        @return: a new pruned catalog
        """

        TIME_DELTA = 60 * 15 # events must be within 15 minutes of each other
        MAG_DELTA = 0.3      # and close enough in magnitude to qualify for
                             # pruning due to spatial proximity

        newCat = self.cat.copy()

        if(time_range is not None):
            st = UTCDateTime(time_range[0]).timestamp
            et = UTCDateTime(time_range[1]).timestamp

            keep_ids = (newCat['EventOrigintim'] >= st) & (newCat['EventOrigintim'] <= et)
            newCat = newCat[keep_ids]
        # end if

        if(None not in [lon, lat, distance_range]):
            n = len(newCat)
            _, _, distances = self.geod.inv(np.ones(n) * lon,
                                            np.ones(n) * lat,
                                            newCat['lon'], newCat['lat'])
            min_dist = degrees2kilometers(distance_range[0]) * 1e3
            max_dist = degrees2kilometers(distance_range[1]) * 1e3
            keep_ids = (distances >= min_dist) & (distances <= max_dist)
            newCat = newCat[keep_ids]
            #print(len(newCat))
        # end if

        if(mag_range is not None):
            keep_ids = (newCat['Mw'] >= mag_range[0]) & (newCat['Mw'] <= mag_range[1])
            newCat = newCat[keep_ids]
        # end if

        if(depth_range is not None):
            keep_ids = (newCat['dep'] >= depth_range[0]) & (newCat['dep'] <= depth_range[1])
            newCat = newCat[keep_ids]
        # end if

        if(min_areal_separation_km > 0):
            n = len(newCat)
            if(n > 0):
                qr = np.ones(n) * self.EARTH_RADIUS_KM
                qt = np.radians(90 - newCat['lat'])
                qp = np.radians(newCat['lon'])
                qxyz = rtp2xyz(qr, qt, qp)
                tree = cKDTree(qxyz)

                id_lists = tree.query_ball_point(qxyz, min_areal_separation_km)

                otimes = np.array(newCat['EventOrigintim'])
                magnitudes = np.array(newCat['Mw'])
                repeated_ids = np.zeros(n, dtype='?')
                for ids in id_lists:
                    prod = np.array(list(product(ids, ids)))
                    prod = prod[~(prod[:, 0] == prod[:, 1])] # drop duplicates

                    # find indices where proximal events are also within TIME_DELTA
                    prod = prod[(np.fabs(otimes[prod[:, 0]] - otimes[prod[:, 1]]) < TIME_DELTA) &
                                (np.fabs(magnitudes[prod[:, 0]] - magnitudes[prod[:, 1]]) < MAG_DELTA)]
                    repeated_ids[prod.flatten()] = True
                # end for

                #pd.set_option('display.max_columns', None)
                #print(newCat[repeated_ids])
                newCat = newCat[~repeated_ids] # drop proximal events
            # end if
        # end if

        return GCMTCatalog(newCat)
    # end func

    def get_origin_timestamps(self):
        return np.array(self.cat['EventOrigintim'])
    # end func

    def get_event_longitudes(self):
        return np.array(self.cat['lon'])
    # end func

    def get_event_latitudes(self):
        return np.array(self.cat['lat'])
    # end func

    def get_event_depths_km(self):
        return np.array(self.cat['dep'])
    # end func

    def get_event_magnitudes(self):
        return np.array(self.cat['Mw'])
    # end func

    def to_obspy_catalog(self):
        return Catalog([event for event in self])
    # end func

    def __iter__(self):
        for i in np.arange(len(self.cat)):
            row = self.cat.iloc[i]

            event = Event(event_type="earthquake",
                          creation_info="")
            origin = Origin()
            magnitude = Magnitude()

            origin.time = row['EventOrigintim']
            origin.latitude = row['lat']
            origin.longitude = row['lon']
            origin.depth = row['dep'] * 1e3 # in m

            magnitude.mag = row['Mw']
            magnitude.magnitude_type = "Mw"

            event.origins.append(origin)
            event.preferred_origin_id = origin.resource_id
            event.magnitudes.append(magnitude)
            event.preferred_magnitude_id = magnitude.resource_id
            yield event
        # end for
    # end func

    def get_compatible_events(self, station_lon1, station_lat1,
                              station_lon2, station_lat2,
                              max_areal_separation_km=15,
                              max_depth_separation_km=10,
                              min_magnitude=-1,
                              az_range=80, min_event_dist_deg=10,
                              max_event_dist_deg=100,
                              max_mt_angle=15):
        """
        For a given pair of stations, finds a list of pairs of proximal earthquakes that meet the
        given criteria for event proximity, distance, azimuth and magnitude
        @param station_lon1: longitude of station 1
        @param station_lat1: latitude of station 1
        @param station_lon2: longitude of station 2
        @param station_lat2: latitude of station 2
        @param max_areal_separation_km: maximum areal separation of event-pair
        @param max_depth_separation_km: maximum separation of event-pair in depth
        @param min_magnitude: minimum magnitude of earthquakes
        @param az_range: (+/-) azimuth range
        @param min_event_dist_deg: minimum distance of event epicentres from either station in degrees
        @param max_event_dist_deg: maximum distance of event epicentres from either station in degrees
        @param max_mt_angle: maximum moment-tensor angle between event-pairs
        @return: dictionary indexed by a pair of event IDs, with the moment-tensor angle between
                 them as the value
        """

        MIN_DIST = degrees2kilometers(min_event_dist_deg) * 1e3
        MAX_DIST = degrees2kilometers(max_event_dist_deg) * 1e3

        p1 = [station_lon1, station_lat1]
        p2 = [station_lon2, station_lat2]
        az, baz, dist = self.geod.inv(p1[0], p1[1], p2[0], p2[1])

        eaz1, ebaz1, edist1 = self.geod.inv(np.ones(len(self.cat)) * p1[0],
                                            np.ones(len(self.cat)) * p1[1],
                                            self.cat['lon'], self.cat['lat'])
        eaz2, ebaz2, edist2 = self.geod.inv(np.ones(len(self.cat)) * p2[0],
                                            np.ones(len(self.cat)) * p2[1],
                                            self.cat['lon'], self.cat['lat'])

        # find event IDs that match the given criteria
        good_ids = ((eaz1 > (baz - az_range)) & (eaz1 < (baz + az_range))) | \
                   ((eaz2 > (az - az_range)) & (eaz2 < (az + az_range)))
        good_ids &= ((edist1 >= MIN_DIST) & (edist1 <= MAX_DIST)) & \
                    ((edist2 >= MIN_DIST) & (edist2 <= MAX_DIST))
        good_ids &= (self.cat['Mw'] >= min_magnitude)

        # find all proximal events
        qr = np.ones(np.sum(good_ids)) * self.EARTH_RADIUS_KM
        qt = np.radians(90 - self.cat['lat'][good_ids])
        qp = np.radians(self.cat['lon'][good_ids])

        qxyz = rtp2xyz(qr, qt, qp)
        id_lists = self.tree.query_ball_point(qxyz, max_areal_separation_km)

        result = defaultdict(list)
        # find angles between all proximal events
        for ids in id_lists:
            if (len(ids) < 2): continue

            ids = np.array(ids)
            # drop events by magnitude
            ids = ids[(self.cat['Mw'][ids] >= min_magnitude)]

            gmat = np.array(self.cat.iloc[ids, 4:10]).T
            gmat_norm = np.linalg.norm(gmat, axis=0)

            gmat /= gmat_norm
            angles = np.arccos(np.clip(np.matmul(gmat.T, gmat), -1, 1))

            mask = np.zeros(angles.shape)
            mask[np.mask_indices(angles.shape[0], np.tril)] = 1
            angles = np.degrees(np.ma.masked_array(angles, mask=mask))

            s_ids = np.argsort(angles.flatten())
            s_ids_i, s_ids_j = np.unravel_index(s_ids, angles.shape)

            for i, j in zip(s_ids_i, s_ids_j):
                if (not np.ma.is_masked(angles[i, j])):
                    # drop events by depth-difference limit
                    depth_difference = np.fabs(self.cat['dep'][ids[i]] - self.cat['dep'][ids[j]])
                    if(depth_difference > max_depth_separation_km): continue

                    if(angles[i, j] > max_mt_angle): break
                    result[(ids[i], ids[j])] = angles[i, j]
                # end if
            # end for
        # end for

        return result
    # end func

    @staticmethod
    def get_mt_angle(e1:Series, e2:Series):
        """
        @param e1: an event row from a GCMTCatalog instance
        @param e2: an event row from a GCMTCatalog instance
        @return: moment-tensor angle between the two events
        """

        mt_angle = np.degrees(np.arccos(np.min([np.dot(e1.iloc[4:10],
                                                       e2.iloc[4:10]) / \
                                               (np.linalg.norm(e1.iloc[4:10]) *
                                                np.linalg.norm(e2.iloc[4:10])), 1.])))
        return mt_angle
    # end func
# end class
