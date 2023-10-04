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
from seismic.misc import split_list

def dump_traces(ds, sn_list, start_date, end_date, length, min_length_sec, output_folder):
    """
    Dump mseed traces from an ASDF file in parallel

    :param ds: ASDF Dataset
    :param sn_list: station list to process
    :param start_date: start date
    :param end_date: end date
    :param length: length of each mseed file
    :param output_folder: output folder
    """

    def make_tag(tr):
        # def make_ASDF_tag(ri, tag):
        data_name = "{net}.{sta}.{loc}.{cha}__{start}__{end}".format(
            net=tr.stats.network,
            sta=tr.stats.station,
            loc=tr.stats.location,
            cha=tr.stats.channel,
            start=tr.stats.starttime.strftime("%Y-%m-%dT%H:%M:%S"),
            end=tr.stats.endtime.strftime("%Y-%m-%dT%H:%M:%S"))
        return data_name
    # end func

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
                    logf.write('Failed to read stream: %s\n' % ('.'.join([net, sta])))
                    continue
                # end try

                if(len(st)==0):
                    current_time += length
                    continue
                # end if

                if(min_length_sec):
                    ost = Stream()

                    for tr in st:
                        if (tr.stats.npts * tr.stats.delta > min_length_sec):
                            ost += tr
                        # end if
                    # end for

                    st = ost
                # end if

                for tr in st:
                    fsName = '{}.MSEED'.format(make_tag(tr))
                    try:
                        tr.write(os.path.join(output_folder, fsName), format="MSEED")
                    except:
                        logf.write('Failed to write stream: %s\n'%(fsName))
                        logf.flush()
                    # end try
                    trCounts += 1
                # end for

                #break
                current_time += length
            # wend
        else:
            # dump all traces as they appear within the asdf file
            logf.write('Exporting mseed files for station: %s\n'%(sn))

            sta = None
            try:
                sta = ds.waveforms[sn]
            except Exception as e:
                print(sn, str(e))
            # end try

            if(sta):
                tags = []
                try:
                    tags = sta.list()
                except Exception as e:
                    print(sn, str(e))
                # end try

                for tag in tags:
                    if ('raw_recording' not in tag): continue

                    st=[]
                    try:
                        st = sta[tag]
                    except Exception as e:
                        print(tag, str(e))
                        continue
                    # end try

                    if (min_length_sec):
                        ost = Stream()

                        for tr in st:
                            if (tr.stats.npts * tr.stats.delta > min_length_sec):
                                ost += tr
                            # end if
                        # end for

                        st = ost
                    # end if

                    for tr in st:
                        if(not isinstance(tr, Trace)): continue

                        fsName = '{}.MSEED'.format(make_tag(tr))
                        try:
                            tr.write(os.path.join(output_folder, fsName),
                                     format="MSEED")
                            trCounts += 1
                            logf.write('Wrote waveform : %s\n' % (fsName))
                        except Exception as e:
                            logf.write('Failed to write trace: %s\n'%(fsName))
                            logf.flush()
                    # end for
                # end for
            # end if
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
@click.option('--min-length-sec', type=int, default=250, help="Minimum length in seconds")
def process(input_asdf, output_folder, start_date, end_date, length, min_length_sec):
    """
    INPUT_ASDF: Path to input ASDF file\n
    OUTPUT_FOLDER: Output folder \n

    The script has two modes of operation: \n
    i. parameters --start-date, --end-date and --length (default 1 day) are specified, which results in
       mseed files, each --length seconds long or shorter (but longer than min_length_sec, if specified),
       being output for the time-range specified \n
    ii. parameters --start-date and --end-date are not provided, which results in the dumping of all
        waveforms (longer than min_length_sec, if specified) as they appear in the ASDF file \n
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

        stations = ds.waveforms.list()
        meta = ds.get_all_coordinates()

        proc_stations = split_list(stations, nproc)

        # output station meta-data
        fn = os.path.join(output_folder, 'stations.txt')
        f = open(fn, 'w+')
        f.write('#Station\t\tLongitude\t\tLatitude\n')
        for sn in stations:
            try:
                f.write('%s\t\t%f\t\t%f\n'%(sn, meta[sn]['longitude'], meta[sn]['latitude']))
            except:
                continue
            # end try
        # end for
        f.close()
    # end if

    # broadcast workload to all procs
    proc_stations = comm.bcast(proc_stations, root=0)

    print (proc_stations[rank])
    dump_traces(ds, proc_stations[rank], start_date, end_date, length, min_length_sec, output_folder)

    del ds
# end func

if (__name__ == '__main__'):
    process()
