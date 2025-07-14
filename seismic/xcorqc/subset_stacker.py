from mpi4py import MPI
from seismic.catalogue.gcmt import GCMTCatalog
from seismic.xcorqc.utils import SpooledMatrix
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp2d
from obspy.geodetics import degrees2kilometers
from seismic.xcorqc.utils import read_subset_stacker_config
import os

DEG2KM = degrees2kilometers(1)

class ExtendedTTInterpolator:
    def __init__(self):
        tt_folder = os.path.dirname(os.path.abspath(__file__)) + '/tt/'
        ptable = loadmat(os.path.join(tt_folder, 'P.mat'))
        stable = loadmat(os.path.join(tt_folder, 'S.mat'))

        self.pio = interp2d(ptable['deg'].squeeze(),
                            ptable['depths'].squeeze(), ptable['tt'].T, kind='cubic')
        self.sio = interp2d(stable['deg'].squeeze(),
                            stable['depths'].squeeze(), stable['tt'].T, kind='cubic')
    # end func

    def get_tt(self, phase, dist_deg, depth_km):
        def get_values(io, x, y):
            if (isinstance(x, np.ndarray) and isinstance(y, np.ndarray)):
                assert len(x) == len(y), 'Length of x and y must be the same'
                rv = np.zeros(len(x))
                for i in np.arange(len(x)):
                    rv[i] = io(x[i], y[i])
                # end for
                return rv
            else:
                rv = io(x, y)
                return rv[0]
            # end if
        # end func

        if (phase == 'P'):
            return get_values(self.pio, dist_deg, depth_km)
        elif (phase == 'S'):
            return get_values(self.sio, dist_deg, depth_km)
        else:
            raise ValueError('TT tables for phase {} not available. Aborting..'. \
                             format(phase))
        # end if
    # end func
# end class

class SubsetStacker():
    def __init__(self):
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        # read configuration
        d = read_subset_stacker_config()
        self.gcmt_catalog_fn = d['CMT_CATALOG_PATH']
        self.V_MINMAX = [d['SW_VMIN'], d['SW_VMAX']]  # km/s
        self.DIST_MINMAX = [d['DIST_MIN'], d['DIST_MAX']]
        self.EMAG_MINMAX = [d['EMAG_MIN'], d['EMAG_MAX']]
        self.AZ_TOL = d['AZ_TOL']
        self.param_dict = d

        self.gc = None
        if(self.rank == 0):
            self.gc = GCMTCatalog(self.gcmt_catalog_fn, ellipse='sphere')
        # end if

        # broadcast gcmt-catalog instance to all ranks
        self.comm.Barrier()
        self.gc = self.comm.bcast(self.gc, root=0)

        # instantiate extended travel-time interpolator
        self.tti = ExtendedTTInterpolator()
    # end func

    def stack(self, spooled_matrix:SpooledMatrix,
              window_start_times:np.ndarray,
              windows_end_times: np.ndarray,
              slon1:float, slat1:float, slon2:float, slat2:float):
        """
        @param spooled_matrix: SpooledMatrix containing cross-correlation entries as rows
        @param window_start_times: start-times corresponding to each row in SpooledMatrix
        @param window_end_times: end-time corresponding to each row in SpooledMatrix
        @param slon1: longitude of site 1
        @param slat1: latitude of site 1
        @param slon2: longitude of site 1
        @param slat2: latitude of site 2
        @return: returns np.ndarrays mean, mean_Xei, mean_Xec, mean_XeiUXec and mean_Xeo
                 as defined in the manuscript
        """

        def circular_select(angles, min_angle, max_angle):
            # Normalize angles to [0, 360)
            angles = np.mod(angles, 360);
            min_angle = np.mod(min_angle, 360);
            max_angle = np.mod(max_angle, 360);

            result = None
            if min_angle <= max_angle:
                # Linear range (e.g., 10-20)
                result = (angles >= min_angle) & (angles <= max_angle)
            else:
                # Circular wraparound range (e.g., 350-10)
                result = (angles >= min_angle) | (angles <= max_angle)
            # end if

            return result
        # end func

        def get_affected_indices(source_eids, pat, swat, swet):
            """
            Finds indices of CC windows affected by P and SW energy
            @param source_eids: source event indices, e.g. for station 1 within/outside azimuth
            @param pat: P arrival time for given station
            @param swat: SW arrival time for given station
            @param swet: SW end time for given station
            @return: returns CC window indices affected by P and/or SW energy
            """
            affected_indices = np.zeros(len(xcst), dtype='?')

            for eid in np.where(source_eids)[0]:
                for i in np.where((pat[eid] >= xcst) & (pat[eid] <= xcet))[0]:
                    affected_indices[i] = 1
                    for j in np.where((swat[eid] >= xcst) & (swat[eid] <= xcet))[0]:
                        affected_indices[j] = 1
                        if (j > i): affected_indices[i:j] = 1

                        for k in np.where((swet[eid] >= xcst) & (swet[eid] <= xcet))[0]:
                            affected_indices[k] = 1
                            if (k > j): affected_indices[j:k] = 1
                        # end for
                    # end for
                # end for
            # end for
            return affected_indices
        # end func

        assert (spooled_matrix.nrows == len(window_start_times) == len(windows_end_times)), \
        'Number of rows in the spooled-matrix should equal the number of window start- and ' \
        'end-times provided'

        xcst = window_start_times
        xcet = windows_end_times

        # relevant events
        reids = (self.gc.cat['EventOrigintim'] >= xcst[0]) & \
                (self.gc.cat['EventOrigintim'] <= xcet[-1]) & \
                (self.gc.cat['MwG'] >= self.EMAG_MINMAX[0]) & \
                (self.gc.cat['MwG'] <= self.EMAG_MINMAX[1])
        cat = self.gc.cat[reids] # subset of events relevant for CC being processed

        # compute distances and azimuths of relevant ecents from the two stations
        p1 = [slon1, slat1]
        p2 = [slon2, slat2]
        az, baz, dist = self.gc.geod.inv(p1[0], p1[1], p2[0], p2[1])

        #print(az, baz)
        eaz1, ebaz1, edist1 = self.gc.geod.inv(np.ones(len(cat)) * p1[0],
                                               np.ones(len(cat)) * p1[1],
                                               cat['lon'], cat['lat'])
        eaz2, ebaz2, edist2 = self.gc.geod.inv(np.ones(len(cat)) * p2[0],
                                               np.ones(len(cat)) * p2[1],
                                               cat['lon'], cat['lat'])

        edistkm1 = edist1 / 1e3
        edistkm2 = edist2 / 1e3

        edistdeg1 = edistkm1 / DEG2KM
        edistdeg2 = edistkm2 / DEG2KM

        # compute P and SW arrival times for relevant events at the two stations
        edepth_km = np.array(cat['dep'])
        otime = np.array(cat['EventOrigintim'])

        ptt1 = self.tti.get_tt('P', edistdeg1, edepth_km)
        ptt2 = self.tti.get_tt('P', edistdeg2, edepth_km)

        pat1 = ptt1 + otime
        pat2 = ptt2 + otime

        swat1 = (edistkm1 / self.V_MINMAX[1]) + otime
        swet1 = (edistkm1 / self.V_MINMAX[0]) + otime

        swat2 = (edistkm2 / self.V_MINMAX[1]) + otime
        swet2 = (edistkm2 / self.V_MINMAX[0]) + otime

        # find event indices that meet distance criteria for stations 1 and 2
        eids1 = (edistdeg1 >= self.DIST_MINMAX[0]) & (edistdeg1 <= self.DIST_MINMAX[1])
        eids2 = (edistdeg2 >= self.DIST_MINMAX[0]) & (edistdeg2 <= self.DIST_MINMAX[1])

        # find event indices within azimuth of stations 1 and 2
        eids1_inside_az = eids1 & circular_select(eaz1, (baz - self.AZ_TOL), (baz + self.AZ_TOL))
        eids2_inside_az = eids2 & circular_select(eaz2, (az - self.AZ_TOL), (az + self.AZ_TOL))
        eids_inside_az = eids1_inside_az | eids2_inside_az

        # find event indices outside azimuth of both stations
        eids_outside_az = ~(eids_inside_az)

        if(True):
            # sanity check
            assert len(set(np.where(eids1_inside_az | eids2_inside_az)[0]).intersection( \
                                    set(np.where(eids_outside_az)[0]))) == 0
        # end if

        # find indices of CCs inside/outside azimuth of relevant events
        ccids1_inside_az = get_affected_indices(eids1_inside_az, pat1, swat1, swet1)
        ccids2_inside_az = get_affected_indices(eids2_inside_az, pat2, swat2, swet2)
        ccids_inside_az = ccids1_inside_az | ccids2_inside_az

        #print('ccs inside azimuth for station 1: {}'.format(np.sum(ccids1_inside_az)))
        #print('ccs inside azimuth for station 2: {}'.format(np.sum(ccids2_inside_az)))

        ccids1_outside_az = get_affected_indices(eids_outside_az, pat1, swat1, swet1)
        ccids2_outside_az = get_affected_indices(eids_outside_az, pat2, swat2, swet2)
        ccids_outside_az = ccids1_outside_az | ccids2_outside_az

        #print('ccs outside azimuth for station 1: {}'.format(np.sum(ccids1_outside_az)))
        #print('ccs outside azimuth for station 2: {}'.format(np.sum(ccids2_outside_az)))

        # aliases to indices as named in the manuscript
        idsXei = ccids_inside_az  # inside az
        idsXec = ~(ccids_inside_az | ccids_outside_az)  # no events
        idsXeiUXec = idsXei | idsXec  # no events outside az range
        idsXeo = ccids_outside_az

        # compute means
        mean = mean_Xei = mean_Xec = mean_XeiUXec = mean_Xeo = None
        mean = np.zeros(spooled_matrix.ncols)
        mean_Xei = np.zeros(spooled_matrix.ncols)
        mean_Xec = np.zeros(spooled_matrix.ncols)
        mean_XeiUXec = np.zeros(spooled_matrix.ncols)
        mean_Xeo = np.zeros(spooled_matrix.ncols)

        for i in np.arange(spooled_matrix.nrows):
            row = spooled_matrix.read_row(i)
            mean += row

            if (idsXei[i]): mean_Xei += row
            if (idsXec[i]): mean_Xec += row
            if (idsXeiUXec[i]): mean_XeiUXec += row
            if (idsXeo[i]): mean_Xeo += row
        # end for

        wc = spooled_matrix.nrows
        wc_Xei = np.sum(idsXei)
        wc_Xec = np.sum(idsXec)
        wc_XeiUXec = np.sum(idsXeiUXec)
        wc_Xeo = np.sum(idsXeo)

        if(wc > 0): mean /= float(wc)
        if(wc_Xei > 0): mean_Xei /= float(wc_Xei)
        if(wc_Xec > 0): mean_Xec /= float(wc_Xec)
        if(wc_XeiUXec > 0): mean_XeiUXec /= float(wc_XeiUXec)
        if(wc_Xeo > 0): mean_Xeo /= float(wc_Xeo)

        """
        np.savez('stack3outputs.npz', xcf=mean,
                 xcf1=mean_Xei, xcf2=mean_Xec, xcf3=mean_XeiUXec,
                 xcf4=mean_Xeo, idsXei=idsXei, idsXec=idsXec,
                 idsXeiUXec=idsXeiUXec, idsXeo=idsXeo)
        """

        return mean, mean_Xei, mean_Xec, mean_XeiUXec, mean_Xeo, \
               wc, wc_Xei, wc_Xec, wc_XeiUXec, wc_Xeo
    # end func
# end class

if __name__=="__main__":
    ss = SubsetStacker()

    for i in np.arange(ss.nproc):
        if(i==ss.rank):
            pass
            #print('Rank: {}\n================'.format(i))
            #print(ss.gc.cat['EventOrigintim'])
        # end if
    # end for
# end if