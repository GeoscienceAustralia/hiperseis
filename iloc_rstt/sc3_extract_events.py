#!/usr/env python
"""
Description:
    This script reads output files produced by 'iloc_phase_ident.py',
    dumps out sc3ml files from the Seiscomp3 database, reads them in and generates
    an ensemble text file containing picks that can be used as input for a
    tomographic inversion

References:

CreationDate:   12/13/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     12/12/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import sys
import click
from obspy import Stream, Trace, UTCDateTime
import subprocess
from tqdm import tqdm
from mpi4py import MPI
from collections import defaultdict 
import numpy as np
from random import shuffle
import os, glob, fnmatch, sys
from obspy.core.event import read_events
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
import MySQLdb
import os

DEVNULL = open(os.devnull, 'wb')

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

def recursive_glob(treeroot, pattern):
    results = []
    for base, dirs, files in os.walk(treeroot):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results
# end func

def runprocess(cmd, get_results=False):
    results = []
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=DEVNULL)
    for line in p.stdout:
        if (get_results): results.append(line.strip())
        #else: print line
    # end for

    p.wait()

    return p.returncode, results
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('data-path', required=True,
                type=str)
@click.argument('scratch-path', required=True,
                type=str)
@click.argument('output-file-stem', required=True,
                type=click.Path(exists=False))
def process(data_path, scratch_path, output_file_stem):
    """
    DATA_PATH: input-folder
    SCRATCH_PATH: scratch-folder
    OUTPUT_FILE_STEM: output file stem
    """

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    # Fetch and distribute event-ids
    eventIds = []
    inventory = defaultdict(list)
    if (rank==0):
        files = recursive_glob(data_path, "*.txt")
        
        for f in files:
            d = np.genfromtxt(f, dtype=[('mstring','S100'),('mfloat','f8')], skip_header=1)
            for item in d:
                if(item[1]>0): eventIds.append(item[0])
            # end for
        # end for
        
        s = set(eventIds)
        assert len(s) == len(eventIds), 'Duplicate event-ids found'
        print 'Processing %d events..'%(len(eventIds)) 
        
        shuffle(eventIds)
        eventIds = split_list(eventIds, nproc)

        # fetch inventory
        db = None
        try:
            db = MySQLdb.connect(host="localhost",
                                 user="sysop",
                                 passwd="sysop",
                                 db="seiscomp3")
        except:
            raise Exception('Failed to connect to database')
        # end try
        
        c = db.cursor()
        c.execute('select n.code, s.code, s.longitude, s.latitude from Network n, Station s where s._parent_oid=n._oid group by n.code, s.code')
        rows = c.fetchall()
        for row in rows:
            inventory['%s.%s'%(row[0], row[1])] = [row[2], row[3]]
        # end for

        db.close()
    #end if
    
    eventIds = comm.scatter(eventIds, root = 0)
    inventory = comm.bcast(inventory, root = 0)
    
    pprocfile = open('%s/pproc.%d.txt'%(scratch_path, rank), 'w+')
    sprocfile = open('%s/sproc.%d.txt'%(scratch_path, rank), 'w+')
    for eid in tqdm(eventIds):
        ofn = '%s/%s.xml'%(scratch_path, eid.split('/')[-1])
        cmd = ['scxmldump -d  mysql://sysop:sysop@localhost/seiscomp3 -E %s -PAMf -o %s'% 
                (eid, ofn)]
        rc, _ = runprocess(cmd, get_results=False)
        if(rc!=0):
            print 'Error exporting event: %s'%(eid)
        else:
            notFound = defaultdict(int)
            cat = read_events(ofn, format='SC3ML')
            
            linesp = []
            liness = []
            for e in cat.events:
                po = e.preferred_origin()

                # retrieve depth; some preferred origins don't have depth values
                poDepth = po.depth
                if(poDepth == None):
                    for o in e.origins:
                        if(o.depth): poDepth = o.depth
                    # end for
                # end if

                if(poDepth == None): continue
                
                for a in po.arrivals:
                    try:
                        ncode = a.pick_id.get_referred_object().waveform_id.network_code
                        scode = a.pick_id.get_referred_object().waveform_id.station_code
                        ccode = a.pick_id.get_referred_object().waveform_id.channel_code
                    except:
                        continue

                    
                    slon = None
                    slat = None
                    try:
                        slon, slat = inventory['%s.%s'%(ncode, scode)]
                    except:
                        notFound['%s.%s'%(ncode, scode)] += 1
                        continue
                    # end try
                    
                    # get band-index and snr from comments
                    pick_attribs = defaultdict(float)
                    pick = a.pick_id.get_referred_object()
                    for c in pick.comments:
                        if ('text' in c.keys()):
                            for item in c['text'].split(','):
                                k,v = item.split('=')
                                pick_attribs[k.strip()] = float(v)
                            # end for
                        # end if
                    # end for

                    # get az, baz and distance
                    da = gps2dist_azimuth(po.latitude, po.longitude, slat, slon)

                    # create row
                    if(a.phase not in ['P', 'S']): continue
                    line = [eid, 
                            po.time.timestamp,
                            e.magnitudes[0].mag if (len(e.magnitudes)) else -999,
                            po.longitude,
                            po.latitude,
                            poDepth/1e3,
                            ncode,
                            scode,
                            ccode,
                            pick.time.timestamp,
                            a.phase,
                            slon,
                            slat,
                            da[1], 
                            da[2], 
                            kilometers2degrees(da[0]/1e3),
                            a.time_residual,
                            pick_attribs['phasepapy_snr'],
                            pick_attribs['quality_measure_cwt'],
                            pick_attribs['dom_freq'],
                            pick_attribs['quality_measure_slope'],
                            int(pick_attribs['band_index']),
                            int(pick_attribs['nsigma'])]

                    if(a.phase == 'P'): linesp.append(line)
                    elif(a.phase == 'S'): liness.append(line)
                # end for
            # end for

            for line in linesp:
                pprocfile.write(' '.join([str(item) for item in line]) + '\n')
            # end for
            
            for line in liness:
                sprocfile.write(' '.join([str(item) for item in line]) + '\n')
            # end for
            if (len(notFound)): print 'Rank: %d'%(rank), notFound
        # end if
        if (os.path.exists(ofn)): os.remove(ofn)
        #break
    # end for
    pprocfile.close()
    sprocfile.close()

    header = '#eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp phase stationLon stationLat az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma\n'
    comm.barrier()
    if(rank == 0):
        ofp = open(output_file_stem + '.p.txt', 'w+')
        ofs = open(output_file_stem + '.s.txt', 'w+')
        
        ofp.write(header)
        ofs.write(header)
        for i in range(nproc):
            pfn = '%s/pproc.%d.txt'%(scratch_path, i)
            sfn = '%s/sproc.%d.txt'%(scratch_path, i)
            
            lines = open(pfn, 'r').readlines()
            for line in lines:
                ofp.write(line)
            # end for
            
            lines = open(sfn, 'r').readlines()
            for line in lines:
                ofs.write(line)
            # end for

            if (os.path.exists(pfn)): os.remove(pfn)
            if (os.path.exists(sfn)): os.remove(sfn)
        # end for
        ofp.close()
        ofs.close()
    # end if
# end func

if __name__=="__main__":
    process()
# end if

