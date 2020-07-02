#!/bin/env python
"""
Description:
    Small utility for exporting mseed files from an asdf file in parallel.

References:

CreationDate:   06/08/18
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     04/12/17   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
import os, sys

from os.path import join, exists
from collections import defaultdict
import numpy as np
from obspy import Stream, Trace, UTCDateTime
from multiprocessing import Pool, TimeoutError
import pyasdf
from obspy.core.trace import Trace
import click

def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
# end func

def dump_traces(ds, sn_list, start_date, end_date, length, output_folder):
    """
    Dump mseed traces from an ASDF file in parallel

    :param ds: ASDF Dataset
    :param sn_list: station list to process
    :param start_date: start date
    :param end_date: end date
    :param length: length of each mseed file
    :param output_folder: output folder
    """

    for sn in sn_list:
        logf = open(os.path.join(output_folder, '%s.log.txt'%(sn)), "w+")

        trCounts = 0
        if(start_date and end_date):
            # dump traces within given time-range and of given length each
            current_time = start_date

            logf.write('Exporting mseed files for station: %s\n'%(sn))
            while(current_time < end_date):

                if(end_date - current_time < length):
                    length = end_date - current_time

                net, sta = sn.split('.')

                st = None
                try:
                    st = ds.get_waveforms(net, sta, '*', '*', current_time, current_time+length, tag='raw_recording')
                except:
                    logf.write('Failed to read stream: %s\n' % (fsName))
                    continue
                # end try

                if(len(st)==0):
                    current_time += length
                    continue
                # end if

                fsName = '%s.%s-%s.MSEED'%(sn, current_time.timestamp, (current_time+length).timestamp)

                try:
                    st.write(os.path.join(output_folder, fsName), format="MSEED")
                    trCounts += len(st)
                except:
                    logf.write('Failed to write stream: %s\n'%(fsName))
                    logf.flush()
                # end try

                #break
                current_time += length
            # wend
        else:
            # dump all traces as they appear within the asdf file
            logf.write('Exporting mseed files for station: %s\n'%(sn))
            sta = ds.waveforms[sn]
            for tag in sta.list():
                s = sta[tag]

                for t in s:
                    if(not isinstance(t, Trace)): continue

                    fsName = '%s.%s-%s.MSEED'%(t.id,
                                               t.stats.starttime.strftime("%y-%m-%d.T%H:%M:%S"),
                                               t.stats.endtime.strftime("%y-%m-%d.T%H:%M:%S"))
                    try:
                        t.write(os.path.join(os.path.join(output_folder, sn), fsName),
                                 format="MSEED")
                        trCounts += len(s)
                        logf.write('Wrote waveform with tag: %s\n' % (tag))
                    except Exception as e:
                        logf.write('Failed to write trace: %s\n'%(fsName))
                        logf.flush()
                # end for
            # end for
        # end if

        logf.write('\t Exported (%d) traces.\n' % (trCounts))

        logf.flush()
        logf.close()
    # end for
# end func


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('input-asdf', required=True,
                type=click.Path(exists=True))
@click.argument('output-folder', required=True,
                type=click.Path(exists=True))
@click.option('--start-date', default=None,
              help="Start date-time in UTC format. If specified, both 'end-date' and 'length' must be specified; " 
                   "otherwise this parameter is ignored.",
              type=str)
@click.option('--end-date', default=None,
              help="End date-time in UTC format. If specified, both 'start-date' and 'length' must be specified; " 
                   "otherwise this parameter is ignored.",
              type=str)
@click.option('--length', default=86400,
              help="Length of each trace in seconds. If specified, both 'start-date' and 'end-date' must be specified; " 
                   "otherwise this parameter is ignored; default is 86400 s (1 day)")
def process(input_asdf, output_folder, start_date, end_date, length):
    """
    INPUT_ASDF: Path to input ASDF file\n
    OUTPUT_FOLDER: Output folder \n

    The script has two modes of operations: \n
    i. parameters --start-date, --end-date and --length (default 1 day) are specified, which results in
       mseed files, each --length seconds long, being output for the time-range specified \n
    ii. parameters --start-date and --end-date are not provided, which results in the dumping of all
        waveforms (of any arbitrary lengths) as they appear in the ASDF file \n
    \n
    Example usage:
    mpirun -np 2 python asdf2mseed.py ./test.asdf /tmp/output --start-date 2013-01-01T00:00:00 --end-date 2016-01-01T00:00:00

    Note:
    """

    try:
        start_date = UTCDateTime(start_date) if start_date else None
        end_date   = UTCDateTime(end_date) if end_date else None
        length     = int(length)
    except:
        assert 0, 'Invalid input'
    # end try

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_stations = None
    ds = pyasdf.ASDFDataSet(input_asdf, mode='r')

    if(rank == 0):
        # split work over stations

        stations = list(ds.get_all_coordinates().keys())
        meta = ds.get_all_coordinates()

        proc_stations = split_list(stations, nproc)

        # output station meta-data
        fn = os.path.join(output_folder, 'stations.txt')
        f = open(fn, 'w+')
        f.write('#Station\t\tLongitude\t\tLatitude\n')
        for sn in stations:
            f.write('%s\t\t%f\t\t%f\n'%(sn, meta[sn]['longitude'], meta[sn]['latitude']))
        # end for
        f.close()
    # end if

    # broadcast workload to all procs
    proc_stations = comm.bcast(proc_stations, root=0)

    print (proc_stations[rank])
    dump_traces(ds, proc_stations[rank], start_date, end_date, length, output_folder)

    del ds
# end func

if (__name__ == '__main__'):
    process()
