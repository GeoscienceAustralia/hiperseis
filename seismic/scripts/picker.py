import click
import logging
from obspy.core import read as obspy_read, Stream

import seismic
from seismic import pslog, config

log = logging.getLogger(__name__)


@click.group()
@click.option('-v', '--verbosity',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Level of logging')
def cli(verbosity):
    seismic.pslog.configure(verbosity)


@cli.command()
@click.argument('config_file')
def pick(config_file):
    log.info('Reading config file...')
    cf = seismic.config.Config(config_file)
    st = Stream()
    for f in cf.miniseeds:
        st += obspy_read(f)
