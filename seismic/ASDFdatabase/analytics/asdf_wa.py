#!/bin/env python
"""
Description:
    Generates a waveform-analytics report on raw data in asdf format

References:

CreationDate:   18/01/2024
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     07/02/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""
import os, sys

from obspy import UTCDateTime
import click
from collections import defaultdict
from seismic.ASDFdatabase.utils import MAX_DATE, MIN_DATE
from seismic.ASDFdatabase.analytics.station_analytics import StationAnalytics
from seismic.ASDFdatabase.analytics.utils import get_response
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
import logging
from mpi4py import MPI

logging.basicConfig()

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('network', required=True,
                type=str)
@click.argument('station', required=True,
                type=str)
@click.argument('channel', required=True,
                type=str)
@click.argument('response-database', required=True,
                type=click.Path(exists=True))
@click.argument('output-folder', required=True,
                type=click.Path(exists=True))
@click.option('--start-date', type=str, default=None, show_default=True,
              help="Start date in UTC format for processing data")
@click.option('--end-date', type=str, default=None, show_default=True,
              help="End date in UTC format for processing data")
def process_asdf(asdf_source, network, station, channel, response_database,
                 output_folder, start_date, end_date):
    """
    ASDF_SOURCE: Path to text file containing paths to ASDF files\n
    NETWORK: network code\n
    STATION: station code\n
    CHANNEL: channel code\n
    RESPONSE_DATABASE: Path to database containing response xmls \n
    OUTPUT_FOLDER: Path to output folder\n
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    try:
        start_date = UTCDateTime(start_date) if start_date else None
        end_date   = UTCDateTime(end_date) if end_date else None
    except Exception as e:
        print(str(e))
        raise RuntimeError('Invalid start- or end-dates')
    # end try

    if(start_date and end_date and ((end_date.date - start_date.date).days<=0)):
        raise RuntimeError('Invalid start- and end-dates. Aborting..')
    # end if

    # instantiate FederatedASDFDataSet
    fds = FederatedASDFDataSet(asdf_source)

    sd = MIN_DATE if start_date is None else start_date
    ed = MAX_DATE if end_date is None else end_date
    meta_list = fds.get_stations(sd, ed, network=network, station=station, channel=channel)

    print(meta_list)
    net, sta, loc, cha = meta_list[0][:4]

    print('Loading response database..')
    resp = get_response(response_database, net, sta, loc, cha)
    if(resp is not None): print('Found response: {}'.format(resp))
    else: raise(RuntimeError('No instrument response found. Aborting..'))

    def get_waveforms_func(net, sta, loc, cha, st, et):
        return fds.get_waveforms(net, sta, loc, cha, st, et, nearest_sample=False)
    # end func

    def get_time_range_func(net, sta, loc, cha):
        return fds.get_recording_timespan(net, sta, loc, cha)
    # end func

    sa = StationAnalytics(get_time_range_func, get_waveforms_func,
                          resp, output_folder,
                          None, sd, ed, nproc=1)

    sa.analyse_data(net, sta, loc, cha)
    # end for

    ofn = os.path.join(output_folder, '{}.{}.{}.{}.pdf'.format(net, sta, loc, cha))
    sa.process_results(ofn, net, sta, loc, cha)
    print('Done..')
# end func

if __name__ == "__main__":
    process_asdf()
# end func
