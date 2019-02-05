from __future__ import print_function

import logging

import click
from obspy.core import read as obspy_read, Stream
from seismic.pickers import pickermaps

import seismic
from seismic import config
from seismic.traveltime import pslog

log = logging.getLogger(__name__)


@click.group()
@click.option('-v', '--verbosity',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Level of logging')
def cli(verbosity):
    seismic.traveltime.pslog.configure(verbosity)


@cli.command()
@click.argument('config_file')
def pick(config_file):
    """
    :param config_file: user supplied config file for picking
    :return: tba
    """
    log.info('Reading config file...')
    cf = config.Config(config_file)

    log.info('Preparing time series')
    st = Stream()

    if cf.seeds:
        for f in cf.miniseeds:
            st += obspy_read(f)
        log.info('Miniseeds accumulated')
    else:
        raise NotImplementedError
    log.info('Applying picking algorithm')
    picker = pickermaps[cf.picker['algorithm']](**cf.picker['params'])

    event = picker.event(st, config=cf)
    event.write(filename='test.xml', format='SC3ML')


@cli.command()
@click.argument('config_file')
def locate(config_file):
    """
    :param config_file: user supplied config file for picking
    :return: tba
    """
    log.info('Reading config file...')
    cf = config.Config(config_file)
