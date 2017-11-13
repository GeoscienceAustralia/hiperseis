"""
Clustering of events and station for 3d inversion input files.
"""
from __future__ import print_function
import os
import click
import logging
import csv
from collections import namedtuple
import pandas as pd
from obspy import read_events
import seismic
from seismic import pslog


log = logging.getLogger(__name__)

Station = namedtuple('Station', 'station_code, latitude, longitude, '
                                'elevation, network_code')


@click.command()
@click.argument('events_dir',
                type=click.Path(exists=True, file_okay=False, dir_okay=True,
                                writable=False, readable=True,
                                resolve_path=True))
@click.argument('station_metadata',
                type=click.File(mode='r'))
@click.option('-o', '--output_file',
              type=click.File(mode='w'), default='cluster_p.csv',
              help='Output P picks file')
@click.option('-x', '--nx', type=int, default=1440,
              help='number of segments from 0 to 360 degrees for longitude')
@click.option('-y', '--ny', type=int, default=720,
              help='number of segments from 0 to 180 degrees for latitude')
@click.option('-z', '--dz', type=float, default=25.0,
              help='unit segment length of depth in meters')
@click.option('-v', '--verbosity',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Level of logging')
def cluster(events_dir, station_metadata, output_file, nx, ny, dz, verbosity):

    seismic.pslog.configure(verbosity)

    events = read_events(os.path.join(events_dir, '*.xml')).events

    stations = _read_stations(station_metadata)

    writer = csv.writer(output_file)

    for e in events:
        process_event(e, stations, writer, nx, ny, dz)


@click.command()
@click.argument('output_file',
                type=click.File(mode='r'))
@click.option('-s', '--sorted_file',
              type=click.File(mode='w'), default='sorted.csv',
              help='output sorted and filter file.')
def sort(output_file, sorted_file):
    """
    Sort and filter based on the source and station block number.

    If there are multiple source and station block combinations, we keep the
    row corresponding to the median observered travel time (observed_tt).

    :param output_file: output file from the cluster stage
    :return: None

    """
    column_names = ['source_block', 'station_block',
                    'residual', 'event_number',
                    'source_longitude', 'source_latitude',
                    'source_depth', 'station_longitude',
                    'station_latitude', 'observed_tt', 'P_or_S']

    cluster_data = pd.read_csv(output_file, header=None,
                               names=column_names)
    cluster_data.sort_values(by=['source_block', 'station_block'],
                             inplace=True)
    cluster_data.to_csv(sorted_file, index=False, header=False)
    groups = cluster_data.groupby(by=['source_block', 'station_block'])
    keep = []
    for _, group in groups:
        med = group['observed_tt'].median()
        keep.append(group[group['observed_tt'] == med])

    final_df = pd.concat(keep)

    final_df.to_csv(sorted_file)


def process_event(event, stations, writer, nx, ny, dz):

    # use timestamp as the event number
    ev_number = int(event.creation_info.creation_time.timestamp * 1e6)
    origin = event.preferred_origin()

    # other event parameters we need
    ev_latitude = origin.latitude
    ev_longitude = origin.longitude
    ev_depth = origin.depth

    dx = 360. / nx
    dy = 180. / ny
    event_block = _find_block(dx, dy, dz, nx, ny,
                              origin.latitude,
                              origin.longitude,
                              z=origin.depth)
    for arr in origin.arrivals:
        sta_code = arr.pick_id.get_referred_object(
        ).waveform_id.station_code
        sta = stations[sta_code]

        # TODO: use station.elevation information
        station_block = _find_block(dx, dy, dz, nx, ny,
                                    float(sta.latitude), float(sta.longitude),
                                    z=0.0)

        # phase_type == 1 if P else 2
        phase_type = 1 if arr.phase == 'P' else 2
        if phase_type == 2:
            assert arr.phase == 'S'

        writer.writerow([
            event_block, station_block, arr.time_residual,
            ev_number, ev_longitude, ev_latitude, ev_depth,
            sta.longitude, sta.latitude,
            (arr.pick_id.get_referred_object().time.timestamp -
             origin.time.timestamp),
            phase_type])


def _find_block(dx, dy, dz, nx, ny, lat, lon, z):
    y = 90. - lat
    x = lon if lon > 0 else lon + 360.0
    i = round(x / dx) + 1
    j = round(y / dy) + 1
    k = round(z / dz) + 1
    station_block = (k - 1) * nx * ny + (j - 1) * nx + i
    return int(station_block)


def _read_stations(csv_file):
    """
    :param csv_file: str
        csv stations file handle passed in by click
    :return: stations_dict: dict
        dict of stations indexed by station_code for quick lookup
    """
    stations_dict = {}
    reader = csv.reader(csv_file)
    next(reader)  # skip header
    for station in map(Station._make, reader):
        stations_dict[station.station_code] = station
    return stations_dict


