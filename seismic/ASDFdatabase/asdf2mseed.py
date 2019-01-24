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

def dump_traces(ds, sn_list, start_date, end_date, length, output_folder):

    for sn in sn_list:
        logf = open(os.path.join(output_folder, '%s.log.txt'%(sn)), "w+")

        if(not os.path.exists(os.path.join(output_folder, sn))):
            os.system('mkdir %s'%os.path.join(output_folder, sn))

        trCounts = np.zeros(3)
        if(start_date and end_date):
            # dump traces within given time-range and of given length each
            current_time = start_date

            logf.write('Exporting mseed files for station: %s\n'%(sn))
            while(current_time < end_date):
                zst = ds.ds_jason_db.fetchDataByTime(ds.ds, sn, ds.zchannel,
                                                     current_time,
                                                     current_time+length)        

                nst = ds.ds_jason_db.fetchDataByTime(ds.ds, sn, ds.nchannel,
                                                     current_time,
                                                     current_time+length)        

                est = ds.ds_jason_db.fetchDataByTime(ds.ds, sn, ds.echannel,
                                                     current_time,
                                                     current_time+length)        

                for i, st in enumerate([zst, nst, est]):
                    if(st is None): continue

                    fsName = '%s.%s.MSEED'%(st.traces[0].id, UTCDateTime(current_time).strftime("%y-%m-%d.T%H:%M:%S"))

                    try:
                        st.merge(method=1, fill_value=0)                
                        st.write(os.path.join(os.path.join(OUTPUT_PATH, sn), fsName), 
                                 format="MSEED")
                        trCounts[i] += 1

                        #st.plot()
                    except:
                        logf.write('Failed to write stream: %s\n'%(fsName))
                        logf.flush()
                # end for
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
                        trCounts[0] += 1
                        logf.write('Wrote waveform with tag: %s\n' % (tag))
                    except Exception as e:
                        logf.write('Failed to write trace: %s\n'%(fsName))
                        logf.flush()
                # end for
            # end for
        # end if
    # end for

    if(start_date and end_date):
        logf.write('\t Exported (%d)Z, (%d)N and (%d)E traces.\n' % (trCounts[0], trCounts[1], trCounts[2]))
    else:
        logf.write('\t Exported %d traces.\n' % (trCounts[0]))
    # end if
    logf.flush()
    logf.close()
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
@click.option('--length', default=3600,
              help="Length of each trace in seconds. If specified, both 'start-date' and 'end-date' must be specified; " 
                   "otherwise this parameter is ignored.")
def process(input_asdf, output_folder, start_date, end_date, length):
    """
    INPUT_ASDF: Path to input ASDF file\n
    OUTPUT_FOLDER: Output folder \n

    Example usage:
    mpirun -np 2 python asdf2mseed.py ./test.asdf /tmp/output --start-date 2013-01-01T00:00:00 --end-date 2016-01-01T00:00:00

    Note: 
    """

    try:
        start_date = UTCDateTime(start_date).timestamp if start_date else None
        end_date   = UTCDateTime(end_date).timestamp if end_date else None
        length      = int(length)
    except:
        assert 0, 'Invalid input'
    # end try

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_stations = defaultdict(list)
    ds = pyasdf.ASDFDataSet(input_asdf)

    if(rank == 0):
        # split work over stations in ds1 for the time being

        stations = ds.get_all_coordinates().keys()
        meta = ds.get_all_coordinates()
        count = 0
        for iproc in np.arange(nproc):
            for istation in np.arange(np.divide(len(stations), nproc)):
                proc_stations[iproc].append(stations[count])
                count += 1
        # end for
        for iproc in np.arange(np.mod(len(stations), nproc)):
            proc_stations[iproc].append(stations[count])
            count += 1
        # end for

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
    
    print proc_stations
    dump_traces(ds, proc_stations[rank], start_date, end_date, length, output_folder)
    # end if
# end func

if (__name__ == '__main__'):
    process()
