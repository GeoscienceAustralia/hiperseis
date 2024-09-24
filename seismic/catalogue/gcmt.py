import numpy as np
import pandas as pd
from pandas.core.series import Series
from scipy.spatial import cKDTree
from pyproj import Geod
from obspy.core import UTCDateTime
from collections import defaultdict
from obspy.geodetics.base import degrees2kilometers
from seismic.misc import rtp2xyz

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
    def __init__(self, fn, ellipse='WGS84'):
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

        self.fn = fn
        self.cat = read_gcmt_catalog(fn)
        self.EARTH_RADIUS_KM = 6371.
        self.tree = None
        self.geod = Geod(ellps=ellipse)
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

        # create kDTree for spatial queries if not done so yet
        if(self.tree is None):
            # create kdtree for spatial queries
            r = np.ones(len(self.cat)) * self.EARTH_RADIUS_KM
            t = np.radians(90 - self.cat['lat'])
            p = np.radians(self.cat['lon'])

            xyz = rtp2xyz(r, t, p)
            self.tree = cKDTree(xyz)
        # end if

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
