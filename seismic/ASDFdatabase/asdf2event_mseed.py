#!/bin/env python
"""
Description:
    Small utility for exporting mseed files related to user-supplied events from an asdf file in parallel.

References:

CreationDate:   18/06/21
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     12/06/21   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
import os, sys

from os.path import join, exists
from collections import defaultdict
import numpy as np
from obspy import Stream, Trace, UTCDateTime
from multiprocessing import Pool, TimeoutError
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from obspy.core.trace import Trace
from obspy import Stream, UTCDateTime, read_events
from obspy.geodetics.base import locations2degrees
from obspy.taup import TauPyModel
import click

def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
# end func

def dump_traces(fds, events_xml, sn_list, start_date, end_date, min_dist, max_dist,
                time_before_p, time_after_p, output_folder):
    """
    :param fds: FederatedASDFDataset
    :param events_xml: events catalogue
    :param sn_list: station list to process
    :param start_date: start date
    :param end_date: end date
    :param min_dist: minimum angular distance from event to station
    :param max_dist: maximum angular distance from event to station
    :param output_folder: output folder
    """

    model = TauPyModel(model="iasp91")
    events = read_events(events_xml)
    meta = fds.unique_coordinates

    for sn in sn_list:
        net, sta = sn.split('.')
        logf = open(os.path.join(output_folder, '%s.log.txt'%(sn)), "w+")

        logf.write('Exporting mseed files for station: %s\n' % (sn))
        export_count = 0
        for ev in events:
            po = ev.preferred_origin()
            dist = locations2degrees(meta[sn][1], meta[sn][0], po.latitude, po.longitude)
            if(dist>=min_dist and dist<=max_dist):

                atimes = model.get_travel_times_geo(po.depth/1000., po.latitude,
                                                    po.longitude, meta[sn][1],
                                                    meta[sn][0],
                                                    phase_list=('P',))

                ptime = po.time + atimes[0].time
                stations = fds.get_stations(ptime-time_before_p, ptime+time_after_p,
                                            network=net, station=sta)
                if(len(stations)):
                    st = Stream()
                    for item in stations:
                        subst = fds.get_waveforms(item[0], item[1], item[2], item[3],
                                                  ptime-time_before_p,
                                                  ptime+time_after_p,)
                        for tr in subst: st.append(tr)
                    # end for

                    fname = '%s.%s.%s.%.4d.%.2d.%.2d.%.2d.%.2d.%.2d.mseed'%(item[0], item[1], item[2], \
                                                                         po.time.year,po.time.month,po.time.day, \
		                                                                 po.time.hour,po.time.minute,po.time.second)
                    st.write(os.path.join(output_folder, fname), format='MSEED')
                    print (fname)
                    export_count += 1
                # end if
            # end for
        # end for

        print('%s: Exported (%d) traces.\n' % (sn, export_count))
        logf.write('\t Exported (%d) traces.\n' % (export_count))
        logf.flush()
        logf.close()
    # end for
# end func


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('input-events', required=True,
                type=click.Path(exists=True))
@click.argument('network-name', required=True,
                type=str)
@click.argument('output-folder', required=True,
                type=click.Path(exists=True))
@click.option('--start-date', default='1900-01-01T00:00:00',
              help="Start date-time in UTC format. Default is 1900-01-01T00:00:00",
              type=str)
@click.option('--end-date', default='2200-01-01T00:00:00',
              help="End date-time in UTC format. Default is 2200-01-01T00:00:00",
              type=str)
@click.option('--min-dist', default=30,
              help="Mininum angular distance between event and station. Default is 30")
@click.option('--max-dist', default=90,
              help="Maximum angular distance between event and station. Default in 90")
@click.option('--time-before-p', default=60,
              help="Trace duration before p arrival in seconds. Default is 60s.")
@click.option('--time-after-p', default=120,
              help="Trace duration after p arrival in seconds. Default is 120s")
def process(asdf_source, input_events, network_name, output_folder, start_date, end_date,
            min_dist, max_dist, time_before_p, time_after_p):
    """
    ASDF_SOURCE: Text file containing a list of paths to ASDF files\n
    INPUT_EVENTS: Path to events catalogue in FDSNStationXML format\n
    NETWORK_NAME: Name of network to process \n

    OUTPUT_FOLDER: Output folder \n
    """

    try:
        start_date = UTCDateTime(start_date) if start_date else None
        end_date   = UTCDateTime(end_date) if end_date else None
    except:
        assert 0, 'Invalid input'
    # end try

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_stations = None
    fds = FederatedASDFDataSet(asdf_source, logger=None)

    if(rank == 0):
        # split work over stations

        stations = list(fds.unique_coordinates.keys())
        stations = [s for s in stations if s.split('.')[0] == network_name] # filter stations
        meta = fds.unique_coordinates

        proc_stations = split_list(stations, nproc)

        # output station meta-data
        fn = os.path.join(output_folder, 'stations.txt')
        f = open(fn, 'w+')
        f.write('#Station\t\tLongitude\t\tLatitude\n')
        for sn in stations:
            f.write('%s\t\t%f\t\t%f\n'%(sn, meta[sn][0], meta[sn][1]))
        # end for
        f.close()
    # end if

    # broadcast workload to all procs
    proc_stations = comm.bcast(proc_stations, root=0)

    dump_traces(fds, input_events, proc_stations[rank], start_date, end_date, min_dist, max_dist,
                time_before_p, time_after_p, output_folder)
# end func

if (__name__ == '__main__'):
    process()
