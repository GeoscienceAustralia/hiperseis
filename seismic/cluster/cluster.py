"""
Clustering of events and station for 3d inversion input files.
"""
import os
import click
import logging
from obspy import read_events, read_inventory
import seismic
from seismic import pslog

log = logging.getLogger(__name__)


@click.group()
@click.option('-v', '--verbosity',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Level of logging')
def cli(verbosity):
    seismic.pslog.configure(verbosity)


@cli.command()
@click.argument('events_dir',
                type=click.Path(exists=True, file_okay=False, dir_okay=True,
                                writable=False, readable=True,
                                resolve_path=True))
@click.argument('station_metadata',
                type=click.File(mode='r'))
@click.option('-x', '--nx', type=int, default=1440,
              help='number of segments from 0 to 360 for longtitude')
@click.option('-y', '--ny', type=int, default=720,
              help='number of segments from 0 to 180 for latitude')
@click.option('-z', '--dz', type=float, default=25.0,
              help='unit segment length of depth in meters')
def cluster(events_dir, station_metadata, nx, ny, dz):
    events = read_events(os.path.join(events_dir, '*.xml')).events
    stations = read_inventory(station_metadata)

    for e in events:

        origin = e.preferred_origin()
        dx = 360./nx
        dy = 180./ny

        y = 90. - origin.latitude

        if origin.longitude > 0.:
            x = origin.longitude
        else:
            x = origin.longitude + 360.
        # endif

        z = origin.depth

        i = round(x/dx) + 1
        j = round(y/dy) + 1
        k = round(z/dz) + 1

        event_block = (k-1) * nx * ny + (j-1) * nx + i

        for arr in origin.arrivals:
            y = 90. - slat

            if slon > 0:
                x = slon
            else:
                x = slon+360.
            # endif

            z = 0.

            i = round(x/dx) + 1
            j = round(y/dy) + 1
            k = round(z/dz) + 1
            station_block = (k-1)*nx*ny+(j-1)*nx+i
