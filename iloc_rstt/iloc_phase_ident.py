#!/usr/env python
"""
Description:
    This script runs iLoc in phase-identification mode, in parallel, and can be used
    to update the underlying Seiscomp3 database

References:

CreationDate:   12/12/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     12/12/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import sys
import click
from obspy import Stream, Trace, UTCDateTime
from tqdm import tqdm
from mpi4py import MPI
from collections import defaultdict 
import numpy as np
from random import shuffle
#import subprocess
import subprocess32 as subprocess
import psutil

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

def kill(proc_pid):
    process = psutil.Process(proc_pid)
    for proc in process.children(recursive=True):
        proc.kill()
    # end for
    process.kill()
# end func                           

def runprocess(cmd, get_results=False):
    results = []
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    try:
        pstdout, pstderr = p.communicate(timeout=600)
        for line in pstdout.splitlines():
            if (get_results): results.append(line.strip())
            #else: print line
        # end for
    except subprocess.TimeoutExpired:
        print (('KILLING PROC: %d'%(p.pid)))
        kill(p.pid)
        p.returncode = 1
    # end try


    #for line in pstderr.splitlines():
    #    print line
    #end for

    return p.returncode, results
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('start-time', required=True,
                type=str)
@click.argument('end-time', required=True,
                type=str)
@click.argument('output-path', required=True,
                type=click.Path(exists=True))
@click.option('--update-db', default=0, help='Update database')
@click.option('--limit-events', default=-1, help='Number of events to process')
def process(start_time, end_time, output_path, update_db, limit_events):
    """
    START_TIME: start time e.g. 1990-04-02T00:00:00 \n
    END_TIME: end time \n
    OUTPUT_PATH: output folder \n
    """

    st = None
    et = None
    try:
        st = UTCDateTime(start_time)
        et = UTCDateTime(end_time)
    except:
        raise Exception('Error in time format')
    # end try

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    # Fetch and distribute event-ids
    eventIds = []
    if (rank==0):
        cmd = ['scevtls -d mysql://sysop:sysop@localhost/seiscomp3 --begin "%s" --end "%s"'% 
                (start_time.replace('T', ' '), end_time.replace('T', ' '))]
        rc, eventIds = runprocess(cmd, get_results=True)
        if(rc==0):
            shuffle(eventIds)
            eventIds = split_list(eventIds, nproc)
        # end if
    # end if

    eventIds = comm.scatter(eventIds, root=0)
    
    convergedCount = 0
    # create output files for each rank
    f = open('%s/proc.%d.txt'%(output_path, rank), 'w+')
    f.write('{0:50}{1}\n'.format('#EventID', 'Converged'))
    
    # apply 'limit_events'
    if(limit_events>0): eventIds = eventIds[:limit_events]
    i = 0
    for eventId in tqdm(eventIds):
        converged = False
        cmd = ['echo "%s update_db=%d do_gridsearch=1 use_RSTT_PnSn=1 use_RSTT_PgLg=1 verbose=0" | iloc sc3db' % 
                (eventId, update_db)]
        
        rc, results = runprocess(cmd, get_results=True)
        if(rc): 
            print (('Event: %s has encountered errors'%(eventId)))
            print (results)
        else:
            for r in results: 
                if ('unexpected number' in r): print (('Event: %s has missing stations [%s]'%(eventId, r)))
                elif('converged: 1' in r): 
                    print (('Event: %s converged [%s]'%(eventId, r)))
                    convergedCount += 1
                    converged = True
                elif('no phases found' in r): print (('Event: %s no phases found'%(eventId)))
            # end for
            print (results)
        # end if
        i+=1
        f.write('{0:50}{1}\n'.format(eventId, int(converged)))
    # end for
    print (('rank: %d #events: %d converged: %d'%(rank, len(eventIds), convergedCount)))
    f.close()
# end func

if __name__=="__main__":
    process()
# end if

