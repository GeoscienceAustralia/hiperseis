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

from utils import EventParser, DistAz

def recursive_glob(treeroot, pattern):
    results = []
    for base, dirs, files in os.walk(treeroot):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results
# end func

def split_list(lst, npartitions):
    result = []
    for i in np.arange(npartitions):
        result.append([])
    # end for
    count = 0
    for iproc in np.arange(npartitions):
        for i in np.arange(np.divide(len(lst), npartitions)):
            result[iproc].append(lst[count])
            count += 1
    # end for
    for iproc in np.arange(np.mod(len(lst), npartitions)):
        result[iproc].append(lst[count])
        count += 1
    # end for

    return result
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
    fds = FederatedASDFDataSet(asdf_source, use_json_db=False, logger=None)
    for nc, sc, start_time, end_time in fds.local_net_sta_list():
        day = 24 * 3600
        dayCount = 0
        curr = start_time
        sw_start = datetime.now()
        traceCount = 0
        while (curr < end_time):
            eventIndices = (np.where((allPOTimestamps >= curr.timestamp) & \
                                     (allPOTimestamps <= (curr + day).timestamp)))[0]

            if(eventIndices.shape[0]>0):
                stations = fds.get_stations(curr, curr + day, network=nc, station=sc)
                stations_zch = [s for s in stations if 'Z' in s[3]]  # only Z channels

                #print 'Found %d events to process on rank %d for day: [%f, %f]' % \
                #      (eventIndices.shape[0], rank, curr.timestamp, (curr + day).timestamp)

                for codes in stations_zch:
                    st = fds.get_waveforms(codes[0], codes[1], codes[2], codes[3],
                                           curr,
                                           curr+day, tag='raw_recording', automerge=True)
                    traceCount += len(st)
                # end for
            # end if
            curr += day
            dayCount += 1
        # wend
        sw_stop = datetime.now()
        totalTime = (sw_stop - sw_start).total_seconds()

        print 'Read %d traces on rank %d for network %s station %s in %f seconds'%\
              (traceCount, rank, nc, sc, totalTime), start_time, end_time
    # end for

    del fds
# end func

if (__name__ == '__main__'):
    process()
# end if