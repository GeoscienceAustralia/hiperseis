#!/bin/env python
"""
Description:
    Harvests picks from ASDF data-sets in parallel
References:

CreationDate:   13/09/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     13/09/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
import glob, os, sys
from os.path import join, exists
from collections import defaultdict

from math import radians, cos, sin, asin, sqrt
import numpy as np
import scipy
from scipy.spatial import cKDTree

import xcorqc
from obspy import Stream, Trace, UTCDateTime
import pyasdf
import json
import fnmatch
import operator
from datetime import datetime
from ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

import click
from random import shuffle
from obspy import UTCDateTime, read_events, read_inventory
from obspy.taup.taup_geo import calc_dist
from obspy.clients.iris import Client as IrisClient
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy.signal.trigger import trigger_onset, z_detect, classic_sta_lta, recursive_sta_lta, ar_pick
from obspy.signal.rotate import rotate_ne_rt
from obspy.core.event import Pick, CreationInfo, WaveformStreamID, ResourceIdentifier, Arrival, Event,\
     Origin, Arrival, OriginQuality, Magnitude, Comment

import pywt
from PhasePApy.phasepapy.phasepicker import fbpicker
from PhasePApy.phasepapy.phasepicker import ktpicker
from PhasePApy.phasepapy.phasepicker import aicdpicker

from utils import EventParser, DistAz, Catalog
import psutil
import gc

def extract_p(taupy_model, picker, event, station_longitude, station_latitude,
              st, win_start=-50, win_end=50, resample_hz=20,
              bp_freqmins = [0.5, 2., 5.],
              bp_freqmaxs = [5, 9.9, 9.95],
              margin=10,
              max_amplitude=1e8):

    po = event.preferred_origin
    if(not po): return None

    atimes = []
    try:
        atimes = taupy_model.get_travel_times_geo(po.depthkm, po.lat,
                                                  po.lon, station_latitude,
                                                  station_longitude,
                                                  phase_list=('P',))
    except:
        return None
    # end try

    if (len(atimes) == 0): return None
    tat = atimes[0].time # theoretical arrival time

    try:
        st = st.slice(po.utctime + tat + win_start, po.utctime + tat + win_end)
        st.resample(resample_hz)
    except:
        return None
    # end try

    if(len(st) == 0): return None

    tr = st[0]
    if(type(tr.data) == np.ndarray):
        if(np.max(tr.data) > max_amplitude): return None

        pickslist = []
        snrlist = []
        residuallist = []

        tr.detrend('linear')
        for i in range(len(bp_freqmins)):
            trc = tr.copy()
            trc.filter('bandpass', freqmin=bp_freqmins[i],
                       freqmax=bp_freqmaxs[i], corners=4,
                       zerophase=True)

            try:
                scnl, picks, polarity, snr, uncert = picker.picks(trc)

                for ipick, pick in enumerate(picks):
                    actualArrival = pick - po.utctime
                    residual = tat - actualArrival

                    if (np.fabs(residual) < margin):
                        pickslist.append(pick)
                        snrlist.append(snr[ipick])
                        residuallist.append(residual)

                        #summary = fbpicker.FBSummary(picker, trc)
                        #summary = aicdpicker.AICDSummary(picker, trc)
                        #outputPath = '/g/data/ha3/rakib/seismic/pst/tests/output'
                        #ofn = '%s/%s.%s_%f_%d.png' % (outputPath, scnl, str(po.utctime), snr[0], i)
                        #summary.plot_picks(show=False, savefn=ofn)
                    # end if
                # end for
            except:
                continue
            # end try
        # end for

        if (len(pickslist)):
            optimal_pick_idx = np.argmax(np.array(snrlist))

            return pickslist[optimal_pick_idx], residuallist[optimal_pick_idx], \
                   snrlist[optimal_pick_idx], optimal_pick_idx
        # end if
    # end if

    return None
# end func

def extract_s(taupy_model, picker, event, station_longitude, station_latitude,
              stn, ste, ba, win_start=-50, win_end=50, resample_hz=20,
              bp_freqmins= [0.05, 2],
              bp_freqmaxs = [0.5, 5],
              margin=20,
              max_amplitude=1e8):

    po = event.preferred_origin
    if(not po): return None

    atimes = []
    try:
        atimes = taupy_model.get_travel_times_geo(po.depthkm, po.lat,
                                                  po.lon, station_latitude,
                                                  station_longitude,
                                                  phase_list=('S',))
    except:
        return None
    # end try

    if (len(atimes) == 0): return None
    tat = atimes[0].time # theoretical arrival time

    tr = None
    try:
        stn = stn.slice(po.utctime + tat + win_start, po.utctime + tat + win_end)
        stn.resample(resample_hz)

        if(ste):
            ste = ste.slice(po.utctime + tat + win_start, po.utctime + tat + win_end)
            ste.resample(resample_hz)
        # end if

        if(ste):
            if(type(stn[0].data) == np.ndarray and type(ste[0].data) == np.ndarray):
                rc, tc = rotate_ne_rt(stn[0].data, ste[0].data, ba)
                tr = Trace(data=tc, header=stn[0].stats)
                #tr = Trace(data=np.sqrt(np.power(rc,2) + np.power(tc,2)), header=stn[0].stats)
            # end if
        else:
            if(type(stn[0].data) == np.ndarray):
                tr = stn[0]
            # end if
        # end if
    except Exception as e:
        return None
    # end try

    if(tr):
        if(np.max(tr.data) > max_amplitude): return None

        pickslist = []
        snrlist = []
        residuallist = []

        tr.detrend('linear')
        for i in range(len(bp_freqmins)):
            trc = tr.copy()
            trc.filter('bandpass', freqmin=bp_freqmins[i],
                       freqmax=bp_freqmaxs[i], corners=4,
                       zerophase=True)

            try:
                scnl, picks, polarity, snr, uncert = picker.picks(trc)

                for ipick, pick in enumerate(picks):
                    actualArrival = pick - po.utctime
                    residual = tat - actualArrival

                    if (np.fabs(residual) < margin):
                        pickslist.append(pick)
                        snrlist.append(snr[ipick])
                        residuallist.append(residual)

                        #summary = fbpicker.FBSummary(picker, trc)
                        #summary = aicdpicker.AICDSummary(picker, trc)
                        #outputPath = '/home/rakib/work/pst/picking/sarr'
                        #outputPath = '/g/data1a/ha3/rakib/seismic/pst/tests/plots/new'
                        #ofn = '%s/%s.%s_%f_%d.s.png' % (outputPath, scnl, str(po.utctime), snr[0], i)
                        #summary.plot_picks(show=False, savefn=ofn)
                    # end if
                # end for
            except:
                continue
            # end try
        # end for

        if (len(pickslist)):
            optimal_pick_idx = np.argmax(np.array(snrlist))

            return pickslist[optimal_pick_idx], residuallist[optimal_pick_idx], \
                   snrlist[optimal_pick_idx], optimal_pick_idx
        # end if
    # end if

    return None
# end func

def getWorkloadEstimate(fds, originTimestamps):
    totalTraceCount = 0
    for nc, sc, start_time, end_time in fds.local_net_sta_list():

        day = 24 * 3600
        curr = start_time
        while (curr < end_time):
            eventIndices = (np.where((originTimestamps >= curr.timestamp) & \
                                     (originTimestamps <= (curr + day).timestamp)))[0]

            if(eventIndices.shape[0]>0): totalTraceCount += 1
            curr += day
        # wend
    # end for
    return totalTraceCount
# end func

def dropBogusTraces(st, sampling_rate_cutoff=5):
    badTraces = [tr for tr in st if tr.stats.sampling_rate < sampling_rate_cutoff]

    for tr in badTraces: st.remove(tr)
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('event-folder', required=True,
                type=click.Path(exists=True))
@click.argument('output-path', required=True,
                type=click.Path(exists=True))
def process(asdf_source, event_folder, output_path):
    """
    ASDF_SOURCE: Text file containing a list of paths to ASDF files
    EVENT_FOLDER: Path to folder containing event files\n
    OUTPUT_PATH: Output folder \n
    """

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_workload = None

    if(rank == 0):
        def outputConfigParameters():
            # output config parameters
            fn = 'pick.%s.cfg' % (UTCDateTime.now().strftime("%y-%m-%d.T%H.%M"))
            fn = os.path.join(output_path, fn)

            f = open(fn, 'w+')
            f.write('Parameters Values:\n\n')
            f.write('%25s\t\t: %s\n' % ('ASDF_SOURCE', asdf_source))
            f.write('%25s\t\t: %s\n' % ('EVENT_FOLDER', event_folder))
            f.write('%25s\t\t: %s\n' % ('OUTPUT_PATH', output_path))
            f.close()
        # end func

        outputConfigParameters()
    # end if

    cat = Catalog(event_folder)
    events = cat.get_events()
    originTimestamps = cat.get_preferred_origin_timestamps()

    # ==========================================
    #picker = fbpicker.FBPicker(t_long=5, freqmin=0.1, mode='std', t_ma=20, nsigma=8, \
    #                           t_up=1, nr_len=5, nr_coeff=2, pol_len=10, pol_coeff=10, uncert_coeff=3)
    picker_p = aicdpicker.AICDPicker(t_ma=5, nsigma=8, t_up=1, nr_len=5,
                                   nr_coeff=2, pol_len=10, pol_coeff=10, uncert_coeff=3)
    picker_s = aicdpicker.AICDPicker(t_ma=15, nsigma=8, t_up=1, nr_len=5,
                                   nr_coeff=2, pol_len=10, pol_coeff=10, uncert_coeff=3)

    taupyModel = TauPyModel(model='iasp91')
    fds = FederatedASDFDataSet(asdf_source, use_json_db=False, logger=None)
    workload = getWorkloadEstimate(fds, originTimestamps)

    header = '#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp stationLon stationLat baz distance ttResidual snr bandIndex\n'
    ofnp = os.path.join(output_path, 'p_arrivals.%d.txt' % (rank))
    ofp = open(ofnp, 'w+')
    ofns = os.path.join(output_path, 's_arrivals.%d.txt' % (rank))
    ofs = open(ofns, 'w+')
    ofp.write(header)
    ofs.write(header)

    totalTraceCount = 0
    for nc, sc, start_time, end_time in fds.local_net_sta_list():
        day = 24 * 3600
        dayCount = 0
        curr = start_time
        traceCountP = 0
        pickCountP = 0
        traceCountS = 0
        pickCountS = 0
        sw_start = datetime.now()
        while (curr < end_time):

            eventIndices = (np.where((originTimestamps >= curr.timestamp) & \
                                     (originTimestamps <= (curr + day).timestamp)))[0]

            if(eventIndices.shape[0]>0):
                totalTraceCount += 1
                stations = fds.get_stations(curr, curr + day, network=nc, station=sc)
                stations_zch = [s for s in stations if 'Z' in s[3]]  # only Z channels
                stations_nch = [s for s in stations if 'N' in s[3] or '1' in s[3]]  # only N channels
                stations_ech = [s for s in stations if 'E' in s[3] or '2' in s[3]]  # only E channels

                for codes in stations_zch:
                    st = fds.get_waveforms(codes[0], codes[1], codes[2], codes[3],
                                           curr,
                                           curr+day, tag='raw_recording', automerge=True)

                    dropBogusTraces(st)
                    if (len(st) == 0): continue

                    slon, slat = codes[4], codes[5]
                    for ei in eventIndices:
                        event = events[ei]
                        po = event.preferred_origin
                        da = DistAz(po.lat, po.lon, slat, slon)
                        mag = None
                        if(event.preferred_magnitude): mag = event.preferred_magnitude.magnitude_value
                        elif(len(po.magnitude_list)): mag = po.magnitude_list[0].magnitude_value
                        if(mag == None): mag = np.NaN

                        result = extract_p(taupyModel, picker_p, event, slon, slat, st)
                        if(result):
                            pick, residual, snr, bi = result

                            line = '%s %f %f %f %f %f ' \
                                   '%s %s %s %f %f %f ' \
                                   '%f %f ' \
                                   '%f %f %d\n' % (event.public_id, po.utctime.timestamp, mag, po.lon, po.lat, po.depthkm,
                                                   codes[0], codes[1], codes[3], pick.timestamp, slon, slat,
                                                   da.getBaz(), da.getDelta(),
                                                   residual, snr, bi)
                            ofp.write(line)
                            pickCountP += 1
                        # end if

                        if (len(stations_nch) == 0 and len(stations_ech) == 0):
                            result = extract_s(taupyModel, picker_s, event, slon, slat, st, None, da.getBaz())
                            if (result):
                                pick, residual, snr, bi = result

                                line = '%s %f %f %f %f %f ' \
                                       '%s %s %s %f %f %f ' \
                                       '%f %f ' \
                                       '%f %f %d\n' % (event.public_id, po.utctime.timestamp, mag, po.lon, po.lat, po.depthkm,
                                                       codes[0], codes[1], codes[3], pick.timestamp, slon, slat,
                                                       da.getBaz(), da.getDelta(),
                                                       residual, snr, bi)
                                ofs.write(line)
                                pickCountS += 1
                            # end if
                        # end if
                    # end for

                    traceCountP += len(st)
                # end for

                if(len(stations_nch)>0 and len(stations_nch) == len(stations_ech)):
                    for codesn, codese in zip(stations_nch, stations_ech):
                        stn = fds.get_waveforms(codesn[0], codesn[1], codesn[2], codesn[3],
                                               curr,
                                               curr + day, tag='raw_recording', automerge=True)
                        ste = fds.get_waveforms(codese[0], codese[1], codese[2], codese[3],
                                               curr,
                                               curr + day, tag='raw_recording', automerge=True)

                        dropBogusTraces(stn)
                        dropBogusTraces(ste)

                        if (len(stn) == 0): continue
                        if (len(ste) == 0): continue

                        slon, slat = codesn[4], codesn[5]

                        for ei in eventIndices:
                            event = events[ei]
                            po = event.preferred_origin
                            da = DistAz(po.lat, po.lon, slat, slon)
                            mag = None
                            if (event.preferred_magnitude):
                                mag = event.preferred_magnitude.magnitude_value
                            elif (len(po.magnitude_list)):
                                mag = po.magnitude_list[0].magnitude_value
                            if (mag == None): mag = np.NaN

                            result = extract_s(taupyModel, picker_s, event, slon, slat, stn, ste, da.getBaz())
                            if (result):
                                pick, residual, snr, bi = result

                                line = '%s %f %f %f %f %f ' \
                                       '%s %s %s %f %f %f ' \
                                       '%f %f ' \
                                       '%f %f %d\n' % (event.public_id, po.utctime.timestamp, mag, po.lon, po.lat, po.depthkm,
                                                       codesn[0], codesn[1], '00T', pick.timestamp, slon, slat,
                                                       da.getBaz(), da.getDelta(),
                                                       residual, snr, bi)
                                ofs.write(line)
                                pickCountS += 1
                            # end if
                        # end for

                        traceCountS += (len(stn) + len(ste))
                    # end for
                # end if
            # end if
            curr += day
            dayCount += 1
        # wend
        sw_stop = datetime.now()
        totalTime = (sw_stop - sw_start).total_seconds()

        gc.collect()
        print '(Rank %d: %5.2f%%, %d/%d) Processed %d traces and found %d p-arrivals and %d s-arrivals for ' \
              'network %s station %s in %f s. Memory usage: %5.2f MB.' % \
              (rank, (float(totalTraceCount) / float(workload) * 100) if workload > 0 else 100, totalTraceCount, workload,
               traceCountP + traceCountS, pickCountP, pickCountS, nc, sc, totalTime,
               round(psutil.Process().memory_info().rss / 1024. / 1024., 2))
    # end for
    ofp.close()
    ofs.close()

    print 'Processing complete on rank %d'%(rank)

    del fds
# end func

if (__name__ == '__main__'):
    process()
# end if
