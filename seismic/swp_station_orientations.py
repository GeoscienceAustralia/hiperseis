#!/usr/bin/env python
# coding=utf-8
"""
    Adaptation of DLOPy:

    PRIMARY ORIENTATION PROGRAM
    ADRIAN. K. DORAN
    GABI LASKE
    VERSION 1.0
    REPO: https://github.com/jbrussell/DLOPy_v1.0
    RELEASED APRIL 2017

Reference:
- Doran,Adrian K. et al.
Ocean-Bottom Seismometer Instrument Orientations via Automated Rayleigh-Wave Arrival-Angle Measurements
Bulletin of the Seismological Society of America(2017),107(2):691

"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import numpy.matlib
import numpy as np
import json
import h5py
from collections import defaultdict
from geographiclib.geodesic import Geodesic

import tqdm.auto as tqdm

import seismic.receiver_fn.rf_util as rf_util
from seismic.stream_io import get_obspyh5_index

from obspy.core import UTCDateTime
from obspy.geodetics import gps2dist_azimuth
from seismic.network_event_dataset import NetworkEventDataset
from scipy.stats import circmean as cmean, circstd as cstd
from scipy.stats import hmean as hm
import scipy.signal as sig
import logging
import click
from mpi4py import MPI

logging.basicConfig()

GRV_FN = os.path.join(os.path.dirname(__file__), 'data/grv.h5')

def checklen(st, hrs):
    # checks to see if there is enough downloaded data to run program
    L=len(st)
    for i in np.arange((L)):
        if (UTCDateTime(st[i].stats.endtime)-UTCDateTime(st[i].stats.starttime))+100 < hrs:
            return True
    if np.var(st[0].data)<1 or np.var(st[1].data)<1 or np.var(st[2].data)<1:
        return True
    return False
# end func

# Overall function to get path-averaged group velocity
def pathvels(lat1, lon1, lat2, lon2, map10, map15, map20, map25, map30, map35, map40):
    #### find nearest value in an array
    def nv(x,v):
        # x is array
        # v is value
        idx = (np.abs(x-v)).argmin()
        return x[idx]
    # end func

    Rearth = 6371.25;
    circE = 2 * np.pi * Rearth;

    # Get distance and azimuth
    p = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)

    minor = float(p['s12']) / 1000.00
    major = circE - minor

    l1 = Geodesic.WGS84.Line(lat1, lon1, p['azi1'])
    l2 = Geodesic.WGS84.Line(lat1, lon1, p['azi1'] - 180)

    deg = 111.17
    D1 = np.zeros([1, 2])
    for i in np.arange((361)):
        b = l1.Position(deg * 1000 * i)
        p2 = np.array([b['lat2'], b['lon2']])

        if i == 0:
            D1[0, 0] = p2[0]
            D1[0, 1] = p2[1]
        else:
            D1 = np.vstack((D1, p2))

        bb = Geodesic.WGS84.Inverse(lat2, lon2, b['lat2'], b['lon2'])
        if bb['s12'] <= deg * 1000.0:
            break

    D2 = np.zeros([1, 2])
    for i in np.arange((361)):
        b = l2.Position(deg * 1000 * i)
        p2 = np.array([b['lat2'], b['lon2']])

        if i == 0:
            D2[0, 0] = p2[0]
            D2[0, 1] = p2[1]
        else:
            D2 = np.vstack((D2, p2))

        bb = Geodesic.WGS84.Inverse(lat2, lon2, b['lat2'], b['lon2'])
        if bb['s12'] <= deg * 1000.0:
            break

    """
    We now have lat and lon points along the major and minor great circles.
    We calcaulte the group velocity of at each point, and then
    find the average velocities. 
    """

    for k in np.arange((len(D1))):
        if D1[k, 1] < 0:
            D1[k, 1] += 360
    for k in np.arange((len(D2))):
        if D2[k, 1] < 0:
            D2[k, 1] += 360

    def Ray(D):
        U1 = np.zeros([len(D), 7])
        for k in np.arange((len(D))):
            # do latitude first

            ## get correct precision
            ## designed to match results of Ma et al codes
            if abs(D[k, 1]) < 10:
                D[k, 1] = round(D[k, 1], 5)
            if abs(D[k, 1]) >= 10 and abs(D[k, 1]) < 100:
                D[k, 1] = round(D[k, 1], 4)
            if abs(D[k, 1]) > 100:
                D[k, 1] = round(D[k, 1], 3)
            if abs(D[k, 0]) < 10:
                D[k, 0] = round(D[k, 0], 5)
            if abs(D[k, 0]) >= 10:
                D[k, 0] = round(D[k, 0], 4)
            #
            # find right index
            q = np.where(map10[:, 1] == nv(map10[:, 1], (D[k, 0])))[0]
            qq = np.where(map10[q, 0] == nv(map10[q, 0], (D[k, 1])))[0]
            idx = q[qq]

            # update path
            U1[k, 0] = map10[idx, 2]
            U1[k, 1] = map15[idx, 2]
            U1[k, 2] = map20[idx, 2]
            U1[k, 3] = map25[idx, 2]
            U1[k, 4] = map30[idx, 2]
            U1[k, 5] = map35[idx, 2]
            U1[k, 6] = map40[idx, 2]
        mhz = np.array([10, 15, 20, 25, 30, 35, 40])
        return np.array((mhz, hm(U1, axis=0))).T

    return Ray(D1), Ray(D2)
# end func

# preprocess segments of data
# taper, zerophase filter, detrend
def sw1proc(T,LPF,HPF,corn=4):
    T.taper(type='hann',max_percentage=0.05)
    T.filter("lowpass",freq=LPF,corners=corn,zerophase=True)
    T.filter("highpass",freq=HPF,corners=corn,zerophase=True)
    T.detrend()

    return T
# end func

def getf(freq,A):
    for i in np.arange((len(A))):
        if A[i,0]==freq:
            v=A[i,1]
    return v
# end func

# Resize arrays to all identical shapes
def resiz(x1,x2,x3):
    a1=len(x1); a2=len(x2); a3=len(x3)
    L=min(np.array([a1,a2,a3]))

    return x1[0:L], x2[0:L], x3[0:L]
# end func

# Define rotation function
#   -rotates horizontal components CW from N by alpha (in degrees)
def rot2d(N,E,alpha):
    a=np.deg2rad(alpha)  # convert from degrees to radians
    r=np.cos(a)*N - np.sin(a)*E
    t=np.sin(a)*N + np.cos(a)*E

    return(r,t)
# end func

# DORAN-LASKE calculation for one freq, one orbit of surface wave
def SW1(TT, Rf, LPF, HPF, daz1, A, nameconv, winlen=10.0):
    # event info
    daz = daz1[0] / 1000.0  # convert to KM
    baz = daz1[1]  # angle from station to event

    Rvel = getf(Rf, A)  # Group velocity at Rf
    R1window = (1.0 / (Rf / 1000.0)) * winlen

    # Process
    T = sw1proc(TT.copy(), LPF, HPF)

    # Window info
    arv = 1.0 / Rvel * daz
    r1 = arv - R1window / 2.0
    r2 = arv + R1window / 2.0

    dt = T[0].stats.starttime
    P = T.slice(starttime=dt + r1, endtime=dt + r2)

    rdat = P[0].data
    rdat2 = P[1].data
    vdat = np.imag(sig.hilbert(P[2].data))

    # Ensure all data vectors are same length
    rdat, rdat2, vdat = resiz(rdat, rdat2, vdat)

    # rotate through and find max cc
    degs = 360 * 4
    ang = np.arange((degs))
    cc = np.zeros((degs));
    cc2 = np.zeros((degs));
    cc3 = np.zeros((degs))
    for k in ang:
        n, e = rot2d(rdat, rdat2, k / 4.0)
        covmat = np.corrcoef(n, vdat)
        cc[k] = covmat[0, 1]
        cstar = np.cov(vdat, n) / np.cov(vdat)
        cc2[k] = cstar[0, 1]
        cstar = np.cov(vdat, e) / np.cov(vdat)
        cc3[k] = cstar[0, 1]

    # Keep angle determined by cstar, but use rating from corrcoef
    #   Formulas in Stachnik paper
    ANG = cc2.argmax();  # CC[j]=cc[ANG]
    # correct for angles above 360
    or_ang = (baz - (360 - ANG / 4.0))

    # ADJUST FOR NAMING CONVENTION
    if nameconv == 3: or_ang += 180

    if or_ang < 0: or_ang += 360
    if or_ang >= 360: or_ang -= 360

    return or_ang, cc[ANG]
# end func

def load_grv():
    """
    Load group velocities
    return: dict of numpy array keyed by frequency (mHz)
    """

    result = {}
    hf = h5py.File(GRV_FN, mode='r')

    for k in hf.keys():
        result[k] = np.array(hf[k][:])
    # end for

    return result
# end func

# Organize channels by coordinate system
def org(T,nameconv):
    # Need to add to this in case other type (eg LHZ) is used
    if nameconv==1:
        bhz=T.select(channel="??Z")
        bh1=T.select(channel="??N")
        bh2=T.select(channel="??E")
    if nameconv==2:
        bhz=T.select(channel="??Z")
        bh1=T.select(channel="??1")
        bh2=T.select(channel="??2")
    if nameconv==3:
        bhz=T.select(channel="??Z")
        bh1=T.select(channel="??2")
        bh2=T.select(channel="??1")

    s=bh1+bh2+bhz

    return s
# end func

def compute_phis(ned, grv_dict, logger=None):
    nameconv = 1
    hrs = 60 * 60 * 4

    R1phi = None
    R1cc = None
    R2phi = None
    R2cc = None

    discarded = 0
    sta = None
    # Loop over stations
    for i, (sta, evts) in enumerate(ned.by_station()): # ned contains 1 station
        nevts = len(evts)

        # Initialize surface wave arrays
        numsurfcalcs = 7
        R1phi = np.zeros([nevts, numsurfcalcs]);
        R1cc = np.zeros([nevts, numsurfcalcs])
        R2phi = np.zeros([nevts, numsurfcalcs]);
        R2cc = np.zeros([nevts, numsurfcalcs]);

        # load group velocity maps
        map10 = grv_dict['10']
        map15 = grv_dict['15']
        map20 = grv_dict['20']
        map25 = grv_dict['25']
        map30 = grv_dict['30']
        map35 = grv_dict['35']
        map40 = grv_dict['40']

        # LOOP OVER ALL EVENTS
        for j, (evid, st) in enumerate(evts.items()):
            st = org(st.copy(), nameconv)

            if checklen(st, hrs):
                discarded += 1
                continue
            # end if

            # remove mean and trend
            st.detrend()
            st.detrend('linear')

            sta_lon = st[0].stats.station_longitude
            sta_lat = st[0].stats.station_latitude
            evt_lon = st[0].stats.event_longitude
            evt_lat = st[0].stats.event_latitude

            # get some additional parameters
            daz1 = gps2dist_azimuth(sta_lat, sta_lon, evt_lat, evt_lon)
            daz2 = np.copy(daz1)
            Rearth = 6371.25 * 1000;
            circE = 2 * np.pi * Rearth;
            daz2[0] = circE - daz2[0];
            daz2[1] = daz2[1] + 180  # major & minor arc calculation
            if daz2[1] >= 360: daz2[1] -= 360

            # SURFACE WAVE CALCULATIONS
            # get path-averaged group velocities
            Ray1, Ray2 = pathvels(sta_lat, sta_lon, evt_lat, evt_lon, map10, map15, map20, map25, map30, map35, map40)
            #
            # FOR EACH FREQUENCY AND ORBIT, calculate arrival angle

            ##    # freq 1 (40 mHz)
            Rf = 40.0;
            HPF = 0.035;
            LPF = 0.045

            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz1, Ray1, nameconv, winlen=20.0)
            R1phi[j, 0] = ANG;
            R1cc[j, 0] = cc

            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz2, Ray2, nameconv, winlen=24.0)
            R2phi[j, 0] = ANG;
            R2cc[j, 0] = cc

            ##    # freq 2 (35 mHz)
            Rf = 35.0;
            HPF = 0.030;
            LPF = 0.040
            #
            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz1, Ray1, nameconv, winlen=17.0)
            R1phi[j, 1] = ANG;
            R1cc[j, 1] = cc

            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz2, Ray2, nameconv, winlen=20.0)
            R2phi[j, 1] = ANG;
            R2cc[j, 1] = cc

            #
            ###    # freq 3 (30 mHz)
            Rf = 30.0;
            HPF = 0.025;
            LPF = 0.035
            #
            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz1, Ray1, nameconv, winlen=14.0)
            R1phi[j, 2] = ANG;
            R1cc[j, 2] = cc

            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz2, Ray2, nameconv, winlen=16.0)
            R2phi[j, 2] = ANG;
            R2cc[j, 2] = cc

            # # # freq 4 (25 mHz)
            Rf = 25.0;
            HPF = 0.020;
            LPF = 0.030
            #
            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz1, Ray1, nameconv, winlen=12.0)
            R1phi[j, 3] = ANG;
            R1cc[j, 3] = cc

            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz2, Ray2, nameconv, winlen=13.0)
            R2phi[j, 3] = ANG;
            R2cc[j, 3] = cc

            ###    # freq 5 (20 mHz)
            Rf = 20.0;
            HPF = 0.015;
            LPF = 0.025
            #
            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz1, Ray1, nameconv, winlen=10.0)
            R1phi[j, 4] = ANG;
            R1cc[j, 4] = cc

            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz2, Ray2, nameconv, winlen=10.0)
            R2phi[j, 4] = ANG;
            R2cc[j, 4] = cc

            ###    # freq 6 (15 mHz)
            Rf = 15.0;
            HPF = 0.020;
            LPF = 0.010
            #
            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz1, Ray1, nameconv, winlen=10.0)
            R1phi[j, 5] = ANG;
            R1cc[j, 5] = cc

            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz2, Ray2, nameconv, winlen=10.0)
            R2phi[j, 5] = ANG;
            R2cc[j, 5] = cc

            ###    # freq 7 (10 mHz)
            Rf = 10.0;
            HPF = 0.005;
            LPF = 0.015
            #
            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz1, Ray1, nameconv, winlen=7.0)
            R1phi[j, 6] = ANG;
            R1cc[j, 6] = cc

            ANG, cc = SW1(st.copy(), Rf, LPF, HPF, daz2, Ray2, nameconv, winlen=7.0)
            R2phi[j, 6] = ANG;
            R2cc[j, 6] = cc

            #break
        # end for
    # end for
    if(logger):
        logger.info("Discarded {}/{} events".format(discarded, len(ned.db_sta[sta])))

    nevents = len(ned.db_sta[sta]) - discarded
    return R1cc, R1phi, R2cc, R2phi, nevents
# end func

# keep eqs above certain cc limit
# also keep which earthquakes were kept
# Different version of C1
def C1_2(phi,cc,clim):
    PHI=np.array([]); C=np.array([]); ix=np.array([])
    for i in np.arange((len(phi))):
        if cc[i]>clim:
            PHI=np.append(PHI,phi[i])
            C=np.append(C,cc[i])
            ix=np.append(ix,i)
    return PHI, C, ix
# end func

# median absolute deviation
def mad(x):
    return np.median(abs(x-np.median(x)))
# end func

# Plotting function
def centerat(phi,m=0):
    phinew=np.copy(phi)
    if len(np.shape(phi))==1:
        for i in np.arange((len(phi))):
            if phi[i]>=m+180:
                phinew[i]-=360
            if phi[i]<=m-180:
                phinew[i]+=360
    else:
        for k in np.arange((np.shape(phi)[1])):
            for i in np.arange((np.shape(phi)[0])):
                if phi[i,k]>=m+180:
                    phinew[i,k]-=360
                if phi[i,k]<=m-180:
                    phinew[i,k]+=360
    return phinew
# end func

# reorganize results
def resort(phi,col2):
    phi2=centerat(phi, m=cmean(phi,high=360))
    t=np.zeros((len(phi2),2))
    t[:,0]=phi2; t[:,1]=col2
    t = t[t[:,0].argsort()]
    return t[:,0], t[:,1]
# end func

# remove outliars
def outlier1(Tphi,ix,lim=5.0):
    devs=abs(Tphi-np.median(Tphi))/mad(Tphi)
    ixs=np.where(devs<lim)
    return Tphi[ixs],ix[ixs]
# end func

# bootstrap mean
def boot1(phi,bootnum):
    m=np.zeros((bootnum)); L=len(phi)
    for i in np.arange((bootnum)):
        a=np.random.choice(phi,size=L,replace=True)
        m[i]=cmean(a,high=360)
    return m
# end func

# Get unique events used in final calculation
def uniqueevents(phis,ccs,n,R1cc,R2cc):
    L=len(phis)/2
    ii=np.zeros((len(n)))
    for i in np.arange((len(n))):
        if n[i]<L:
            ii[i]=np.where(R1cc==ccs[int(n[i])])[0][0]
        else:
            ii[i]=np.where(R2cc==ccs[int(n[i])])[0][0]

    return np.unique(ii)
# end func

# final Doran-Laske calculation
def fcalc1(phi,cc,lim,R1cc,R2cc):
    # keep cc over limit
    Tphi,Tcc,ii=C1_2(phi,cc,lim)
    if len(Tphi)==0:
        return 0,180,np.array([0]),0
    if mad(Tphi)==0:
        return np.mean(Tphi),90, np.array([1]),len(Tphi)
    # remove outliers using MAD
    Tphi,ii=resort(Tphi,ii)
    Tphi2,ii2=outlier1(Tphi,ii)
    # bootstrap results for statistic
    m=boot1(Tphi2,5000)


    return cmean(m,high=360),2*1.96*cstd(m,high=360), uniqueevents(phi,cc,ii2,R1cc,R2cc),len(Tphi2)
# end func


def summary_calculations(R1cc, R1phi, R2cc, R2phi, logger=None):
    """
    Final Orientation Calculation File
    A. Doran and G. Laske
    """

    if(logger): logger.info('Calculating summary')

    #####################
    ## ANGLE CALCULATION PARAMETERS
    #####################

    LIM = 0.8  # CC limit for Surface wave calculations

    #
    ### Specify phases to use
    R1use = 1
    R1_40 = 1;
    R1_35 = 1;
    R1_30 = 1;
    R1_25 = 1;
    R1_20 = 1;
    R1_15 = 1;
    R1_10 = 1
    R2use = 1
    R2_40 = 1;
    R2_35 = 1;
    R2_30 = 1;
    R2_25 = 1;
    R2_20 = 1;
    R2_15 = 1;
    R2_10 = 1
    #
    #

    ######################
    ### FINAL ANGLE CALCULATIONS
    ######################
    #
    # Initialize arrays
    L = len(R1phi)
    phis = np.array([])
    ccs = np.array([])
    finval = np.array([]);
    finerr = np.array([])
    N = np.array([]);

    N = np.full((L, L), -1.0)
    LN = np.zeros((L))
    phases = np.array([]);

    startL = 0
    endL = 0

    A = np.array([L])
    if endL != 0:
        A = np.array([endL])

    # If not all calculations are desired, adjust accordingly
    sha = np.shape(R1phi)
    if R1use == 0: R1cc = np.zeros(sha)
    if R1use == 1 and R1_40 == 0: R1cc[:, 0] = np.zeros((sha[0]))
    if R1use == 1 and R1_35 == 0: R1cc[:, 1] = np.zeros((sha[0]))
    if R1use == 1 and R1_30 == 0: R1cc[:, 2] = np.zeros((sha[0]))
    if R1use == 1 and R1_25 == 0: R1cc[:, 3] = np.zeros((sha[0]))
    if R1use == 1 and R1_20 == 0: R1cc[:, 4] = np.zeros((sha[0]))
    if R1use == 1 and R1_15 == 0: R1cc[:, 5] = np.zeros((sha[0]))
    if R1use == 1 and R1_10 == 0: R1cc[:, 6] = np.zeros((sha[0]))

    if R2use == 0: R2cc = np.zeros(sha)
    if R2use == 1 and R2_40 == 0: R2cc[:, 0] = np.zeros((sha[0]))
    if R2use == 1 and R2_35 == 0: R2cc[:, 1] = np.zeros((sha[0]))
    if R2use == 1 and R2_30 == 0: R2cc[:, 2] = np.zeros((sha[0]))
    if R2use == 1 and R2_25 == 0: R2cc[:, 3] = np.zeros((sha[0]))
    if R2use == 1 and R2_20 == 0: R2cc[:, 4] = np.zeros((sha[0]))
    if R2use == 1 and R2_15 == 0: R2cc[:, 5] = np.zeros((sha[0]))
    if R2use == 1 and R2_10 == 0: R2cc[:, 6] = np.zeros((sha[0]))

    # Function to flatten result arrays
    def flatten(X):
        return np.reshape(X,[X.shape[0]*X.shape[1],1])
    # end func

    for i in A:
        # create one massive list with necessary angles and cc values
        phis = np.concatenate((flatten(R1phi[startL:i, :]), flatten(R2phi[startL:i, :])))
        ccs = np.concatenate((flatten(R1cc[startL:i, :]), flatten(R2cc[startL:i, :])))

        # Doran-Laske calculation
        val, err, n, ph = fcalc1(phis, ccs, LIM, R1cc, R2cc)
        finval = np.append(finval, val)
        finerr = np.append(finerr, err)
        phases = np.append(phases, ph)
        for k in np.arange((len(n))):
            N[k, i - 1] = n[k]
        LN[i - 1] = len(n)
    # end for

    #D-L mean, error, data included, unique events
    return finval[-1], finerr[-1], phases[-1], max(LN)
# end func

def analyze_station_orientations(ned, grv_dict, save_plots_path=None, ax=None):
    assert isinstance(ned, NetworkEventDataset), 'Pass NetworkEventDataset as input'

    if len(ned) == 0:
        return {}

    # Determine limiting date range per station
    results = defaultdict(dict)
    full_code = None
    for sta, db_evid in ned.by_station():
        start_time = end_time = None
        full_code = '.'.join([ned.network, sta])
        for _evid, stream in db_evid.items():
            if start_time is None:
                start_time = stream[0].stats.starttime
            # end if
            if end_time is None:
                end_time = stream[0].stats.endtime
            # end if
            for tr in stream:
                start_time = min(start_time, tr.stats.starttime)
                end_time = max(end_time, tr.stats.endtime)
            # end for
        # end for
        results[full_code]['date_range'] = [str(start_time), str(end_time)]
    # end for

    logger = logging.getLogger(__name__ + ':' + full_code)
    logger.setLevel(logging.INFO)

    logger.info('Analysing arrivals')

    r1cc, r1phi, r2cc, r2phi, nevents = compute_phis(ned, grv_dict, logger)
    corr, err, ndata, nevents_c = summary_calculations(r1cc, r1phi, r2cc, r2phi, logger)

    corr *= -1 # converting to azimuth correction
    while (corr > 180): corr -= 360
    while (corr < -180): corr += 360

    results[full_code]['azimuth_correction'] = corr
    results[full_code]['uncertainty'] = err

    logger.info('corr {:2.3f}°, stddev {:2.3f}° (Nc = {:3d})'.format(-corr, err, int(nevents_c)))

    CEN = corr
    c = np.matlib.repmat([40, 35, 30, 25, 20, 15, 10], r1cc.shape[0], 1)
    if save_plots_path is not None:
        _f = plt.figure(figsize=(16, 9))
        plt.subplot(1, 1, 1)
        plt.title('DLOPy results: ' + full_code, fontsize=16)
        plt.plot([0, 1], [corr, corr], '-', linewidth=4, color=(0.8, 0.8, 0.8), zorder=5)
        sc = plt.scatter(r1cc, centerat(-r1phi, m=CEN), c=c, marker='o', cmap=cm.viridis, alpha=0.5, zorder=1,
                         label='R1')
        plt.scatter(r2cc, centerat(-r2phi, m=CEN), c=c, marker='^', cmap=cm.viridis, alpha=0.5, zorder=1,
                    label='R2')
        plt.text(0.5, corr, 'Mean: {0:.3f} deg'.format(corr), fontsize=14, zorder=10)
        cbar = plt.colorbar(sc)
        cbar.set_label('Frequency (mHz)')
        plt.xlabel('Correlation coefficient', fontsize=12)
        plt.ylabel('Orientation correction (deg)', fontsize=12)
        plt.ylim([CEN - 180, CEN + 180]);
        plt.xlim([0, 1])
        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)

        plt.title(full_code, fontsize=14)

        plt.text(0.9, 0.9, 'N = {}'.format(nevents), ha='right', va='top',
                 transform=plt.gca().transAxes)
        plt.text(0.9, 0.8, 'N$_c$ = {}'.format(int(nevents_c)), ha='right', va='top',
                 transform=plt.gca().transAxes)

        plt.legend(framealpha=0.5)
        outfile = '_'.join([full_code, 'swp_ori.png'])
        outfile = os.path.join(str(save_plots_path), outfile)
        plt.savefig(outfile, dpi=300)
        plt.close()
    elif(ax is not None):
        fig = ax.get_figure()
        ax.plot([0, 1], [corr, corr], '-', linewidth=4, color=(0.8, 0.8, 0.8), zorder=5)
        sc = ax.scatter(r1cc, centerat(-r1phi, m=CEN), c=c, marker='o', cmap=cm.viridis, alpha=0.5, zorder=1,
                         label='R1')
        ax.scatter(r2cc, centerat(-r2phi, m=CEN), c=c, marker='^', cmap=cm.viridis, alpha=0.5, zorder=1,
                    label='R2')
        ax.text(0.5, corr, 'Mean: {0:.3f} deg'.format(corr), fontsize=14, zorder=10)

        cax = fig.add_axes([0.1, 0.075, 0.01, 0.1])
        cbar = fig.colorbar(sc, cax=cax, orientation='vertical')
        cbar.set_label('Frequency (mHz)', fontsize=6)
        ax.set_xlabel('Correlation coefficient', fontsize=12)
        ax.set_ylabel('Orientation correction (deg)', fontsize=12)
        ax.set_ylim([CEN - 180, CEN + 180]);
        ax.set_xlim([0, 1])

        ax.text(0.9, 0.9, 'N = {}'.format(nevents), ha='right', va='top',
                transform=ax.transAxes)
        ax.text(0.9, 0.8, 'N$_c$ = {}'.format(int(nevents_c)), ha='right', va='top',
                transform=ax.transAxes)
        ax.legend(framealpha=0.5)
    # end if

    return results
# end func
