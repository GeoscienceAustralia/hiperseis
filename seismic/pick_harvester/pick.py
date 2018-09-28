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

from utils import EventParser, DistAz, recursive_glob, split_list

def extract_p(taupy_model, picker, event, station_longitude, station_latitude,
              st, win_start=-50, win_end=50, resample_hz=20):
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
        st = st.slice(po.utctime + tat + win_start*1.25, po.utctime + tat + win_end*1.25)
        st.resample(resample_hz)
        st = st.slice(po.utctime + tat + win_start, po.utctime + tat + win_end)
    except:
        return None
    # end try

    if(len(st) == 0): return None

    tr = st[0]
    if(type(tr.data) == np.ndarray):
        if(np.max(tr.data) > 1e8): return None

        bp_freqmins = [0.5, 2., 5.]
        bp_freqmaxs = [5, 9.9, 9.95]
        margin = 10
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
                   snrlist[optimal_pick_idx]
        # end if
    # end if

    return None
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
        # retrieve list of all event xml files
        xml_files = recursive_glob(event_folder, '*.xml')
        shuffle(xml_files) # shuffle to avoid concentration of large files
        proc_workload = split_list(xml_files, nproc)
    # end if

    # broadcast workload to all procs
    proc_workload = comm.bcast(proc_workload, root=0)

    print 'Rank %d: processing %d files'%(rank, len(proc_workload[rank]))

    eventList = []
    poTimestamps =[]
    for ifn, fn in enumerate(proc_workload[rank]):
        es = EventParser(fn).getEvents()
        for i, (eid, e) in enumerate(es.iteritems()):
            po = e.preferred_origin
            if(not (po.depthkm >= 0)): continue

            eventList.append(e)
            poTimestamps.append(po.utctime.timestamp)
        # end for
    # end for

    eventList = comm.allgather(eventList)
    poTimestamps = comm.allgather(poTimestamps)

    allEventList = []
    allPOTimestamps = []

    for iproc in np.arange(nproc):
        for i, e in enumerate(eventList[iproc]):
            allEventList.append(e)
            allPOTimestamps.append(poTimestamps[iproc][i])
        # end for
    # end for
    allPOTimestamps = np.array(allPOTimestamps)

    if(rank==0):
        print 'Collected %d event origins'%(len(allEventList))

        hasPM = 0
        hasMultipleMags = 0
        for e in allEventList:
            o = e.preferred_origin
            if(e.preferred_magnitude): hasPM += 1
            if(len(o.magnitude_list)): hasMultipleMags += 1
        # end for

        print '%d preferred origins have a preferred magnitude'%(hasPM)
        print '%d preferred origins have at least one magnitude'%(hasMultipleMags)
    # end if

    # ==========================================
    #picker = fbpicker.FBPicker(t_long=5, freqmin=0.1, mode='std', t_ma=20, nsigma=8, \
    #                           t_up=1, nr_len=5, nr_coeff=2, pol_len=10, pol_coeff=10, uncert_coeff=3)
    picker = aicdpicker.AICDPicker(t_ma=5, nsigma=8, t_up=1, nr_len=5,
                                   nr_coeff=2, pol_len=10, pol_coeff=10, uncert_coeff=3)

    taupyModel = TauPyModel(model='iasp91')
    fds = FederatedASDFDataSet(asdf_source, use_json_db=False, logger=None)
    for nc, sc, start_time, end_time in fds.local_net_sta_list():
        ofn = os.path.join(output_path, '%s.%s.%d.txt'%(nc, sc, rank))
        of = open(ofn, 'w+')

        day = 24 * 3600
        dayCount = 0
        curr = start_time
        sw_start = datetime.now()
        traceCount = 0
        pickCount = 0
        while (curr < end_time):
            eventIndices = (np.where((allPOTimestamps >= curr.timestamp) & \
                                     (allPOTimestamps <= (curr + day).timestamp)))[0]

            if(eventIndices.shape[0]>0):
                stations = fds.get_stations(curr, curr + day, network=nc, station=sc)
                stations_zch = [s for s in stations if 'Z' in s[3]]  # only Z channels

                for codes in stations_zch:

                    #if (codes[0] != 'AU'): continue
                    #if (codes[1] != 'PTPS' or codes[1] != 'AS07'): continue

                    st = fds.get_waveforms(codes[0], codes[1], codes[2], codes[3],
                                           curr,
                                           curr+day, tag='raw_recording', automerge=True)

                    if (len(st) == 0): continue # no data found or likely merge failed
                    if (st[0].stats.sampling_rate < 1): continue # likely corrupt data

                    slon, slat = codes[4], codes[5]
                    for ei in eventIndices:
                        event = allEventList[ei]

                        result = extract_p(taupyModel, picker, event, slon, slat, st)
                        if(result):
                            po = event.preferred_origin
                            da = DistAz(po.lat, po.lon, slat, slon)

                            pick, residual, snr = result
                            #print residual, snr
                            line = '%s %s %s %s %s %f %f %f\n' % (event.public_id, codes[0],
                                                                  codes[1], codes[3], str(pick.time),
                                                                  da.getBaz(), residual, snr)
                            of.write(line)
                            pickCount += 1
                        # end if
                    # end for

                    traceCount += len(st)
                # end for
            # end if
            curr += day
            dayCount += 1
        # wend
        sw_stop = datetime.now()
        totalTime = (sw_stop - sw_start).total_seconds()

        print 'Read %d traces and found %d picks on rank %d for network %s station %s in %f seconds'%\
              (traceCount, pickCount, rank, nc, sc, totalTime), start_time, end_time
        of.close()
    # end for

    del fds
# end func

if (__name__ == '__main__'):
    process()
# end if