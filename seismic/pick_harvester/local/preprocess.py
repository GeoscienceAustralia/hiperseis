#!/usr/bin/env python

"""
Description: This script is largely an adaptation of scmssort, provided by SeisComp3.
Changes were made to make the original scmssort more flexible and fit for our purposes.
The modified script reads in mseed files from an input folder and dumps a single 
multiplexed mseed file. Note that mseed data are sorted in memory and memory requirements
can be astronomical.
    

CreationDate:   03/02/20
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     03/02/20   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import sys, os
from obspy.core import UTCDateTime
import click
from glob import glob
import fnmatch
import numpy as np
from collections import defaultdict
import subprocess32 as subprocess
import seiscomp3.Seismology, seiscomp3.Core, seiscomp3.IO
from tqdm import tqdm 

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
        sys.stderr.write('KILLING PROC: %d'%(p.pid))
        kill(p.pid)
        p.returncode = 1
    # end try

    return p.returncode, results
# end func

def str2time(timestring):
    """
    Liberally accept many time string formats and convert them to a
    seiscomp3.Core.Time
    """
    
    timestring = timestring.strip()
    for c in ["-","/",":", "T", "Z"]:
        timestring = timestring.replace(c, " ")
    timestring = timestring.split()
    assert 3 <= len(timestring) <= 6
    timestring.extend((6-len(timestring))*["0"])
    timestring = " ".join(timestring)
    format = "%Y %m %d %H %M %S"
    if timestring.find(".") != -1:
        format += ".%f"

    t = seiscomp3.Core.Time()
    t.fromString(timestring, format)
    return t
# end func

def _time(rec, sort_by_end_time=True):
    if sort_by_end_time:
        return seiscomp3.Core.Time(rec.endTime())
    return seiscomp3.Core.Time(rec.startTime())
# end func

def _valid_record(rec):
    return rec is not None # may get more complicated ;)
# end func

def RecordInput(filename=None, datatype=seiscomp3.Core.Array.INT):
    """
    Simple Record iterator that reads from a file (to be specified by
    filename) or -- if no filename was specified -- reads from standard input
    """

    stream = seiscomp3.IO.RecordStream.Create("file")
    if not stream:
        raise IOError, "failed 1"

    if not filename:
        filename = "-"

    if not stream.setSource(filename):
        raise IOError, "failed 2"

    input = seiscomp3.IO.RecordInput(
                    stream, datatype, seiscomp3.Core.Record.SAVE_RAW)

    while 1:
        rec = input.next()
        if not rec:
            raise StopIteration

        yield rec
# end func



def get_mseed_files(folder):
    flist = []
    for root, dirnames, filenames in os.walk(folder):
        for filename in fnmatch.filter(filenames, '*MSEED'):
            flist.append(os.path.join(root, filename))
        # end for
    # end for

    stlist = []
    etlist = []
    for f in flist:
        toks = f.split('-')
        st = float('.'.join(toks[0].split('.')[-2:]))
        et = float('.'.join(toks[1].split('.')[0:2]))

        stlist.append(st)
        etlist.append(et)
    # end for
    return np.array(flist), np.array(stlist), np.array(etlist)
# end func    

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('data-folder', 
                type=click.Path('r'))
@click.option('--start-time', default=None, help="Date and time (in UTC format) to start accepting mseed files from")
@click.option('--end-time', default=None, help="Date and time (in UTC format) to stop accepting mseed files until")
@click.option('--span-days', default=7, help="Timespan in days over which to combine mseeds for output; default 7.")
def main(data_folder, start_time, end_time, span_days):
    """
    DATA_FOLDER: Data folder that is recursively scanned for mseed files
    """

    flist, stlist, etlist = get_mseed_files(data_folder)
    
    if(start_time == None and end_time == None):
        st, et = np.min(stlist), np.max(etlist)
    else:
        st = UTCDateTime(start_time).timestamp
        et = UTCDateTime(end_time).timestamp
    # end if
    
    ct = st
    step = span_days*24*3600
    bins = defaultdict(list)
    
    while ct < et:
        if(ct+step > et): step = et-ct

        bins[ct] = flist[(stlist>=ct) & (etlist<=(ct + step))]

        ct += step
    # wend
    
    #assert len(flist) == np.sum([len(v) for k,v in bins.iteritems()]), 'Total # of files do not match total # of binned files..'
    
    step = span_days*24*3600
    for i, (ws,files) in enumerate(bins.iteritems()):
        sys.stderr.write('[%s - %s]: %d files\n'%(UTCDateTime(ws).strftime('%Y-%m-%d'),
                                                UTCDateTime(ws+step).strftime('%Y-%m-%d'),
                                                len(files)))

        #if(i!=5): continue
        
        first = None
        time_raw = []
        for filename in tqdm(files, desc='Reading files', leave=True):
            for rec in RecordInput(filename):
                if not _valid_record(rec):
                    continue

                raw = rec.raw().str()
                t = _time(rec)
                if first is None:
                    first = t
                t = float(t-first) # float needs less memory
                time_raw.append( (t,raw) )
            # end for
        # end for

        time_raw.sort()
        
        ofn = '%s-%s.mseed'%(UTCDateTime(start_time).strftime('%Y.%m.%d'), 
                             UTCDateTime(end_time).strftime('%Y.%m.%d'))
        of = open(ofn, 'wb')
        previous = None
        for item in tqdm(time_raw, desc='Writing records', leave=True):
            if item == previous:
                continue
            t,raw = item        
            of.write(raw)
            previous = item
        # end for
        of.close()
    # end for
# end func


if __name__=='__main__':
    main()
# end if    
