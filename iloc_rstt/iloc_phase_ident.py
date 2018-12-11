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
import subprocess
from tqdm import tqdm
from mpi4py import MPI
from collections import defaultdict 
import numpy as np

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

def runprocess(cmd, get_results=False):
    results = []
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in p.stdout:
        if (get_results): results.append(line.strip())
        #else: print line
    # end for

    #for line in p.stderr:
    #    print line
    p.wait()

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
    START_TIME: start time e.g. 1990-04-02T00:00:00
    END_TIME: end time
    OUTPUT_PATH: output folder
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

    # Fetch event-ids
    cmd = ['scevtls -d mysql://sysop:sysop@localhost/seiscomp3 --begin "%s" --end "%s"'% 
            (start_time.replace('T', ' '), end_time.replace('T', ' '))]
    rc, eventIds = runprocess(cmd, get_results=True)
    
    convergedCount = 0
    if(rc==0):
        # create output files for each rank
        f = open('%s/proc.%d.txt'%(output_path, rank), 'w+')
        f.write('{0:50}{1}\n'.format('#EventID', 'Converged'))
        
        # create sublist of events for each rank and apply 'limit_events'
        eventIds = split_list(eventIds, nproc)[rank]
        if(limit_events>0): eventIds = eventIds[:limit_events]
        i = 0
        for eventId in tqdm(eventIds):
            converged = False
            cmd = ['echo "%s update_db=%d do_gridsearch=0 use_RSTT_PnSn=1 use_RSTT_PgLg=1 verbose=1" | iloc sc3db' % 
                    (eventId, update_db)]
            
            rc, results = runprocess(cmd, get_results=True)
            if(rc): 
                print 'Event: %s has encountered errors'%(eventId)
            else:
                for r in results: 
                    if ('unexpected number' in r): print 'Event: %s has missing stations [%s]'%(eventId, r)
                    elif('converged: 1' in r): 
                        print 'Event: %s converged [%s]'%(eventId, r)
                        convergedCount += 1
                        converged = True
                    elif('no phases found' in r): print 'Event: %s no phases found'%(eventId)
                # end for
            # end if
            i+=1
            f.write('{0:50}{1}\n'.format(eventId, int(converged)))
        # end for
        print 'rank: %d #events: %d converged: %d'%(rank, len(eventIds), convergedCount)
        f.close()
    # end if
# end func

if __name__=="__main__":
    process()
# end if

