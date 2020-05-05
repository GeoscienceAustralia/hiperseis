import logging
import yaml
from os import path
from obspy import UTCDateTime
log = logging.getLogger(__name__)


class Config:
    """Class representing the global configuration of the picking and
    location model to simulate.

    Parameters
    ----------
    yaml_file : string
        The path to the yaml config file. For details on the yaml schema
        see the passive-seismic documentation
    """

    def __init__(self, yaml_file):
        with open(yaml_file, 'r') as f:
            s = yaml.load(f, Loader=yaml.SafeLoader)

        self.name = path.basename(yaml_file).rsplit(".", 1)[0]
        self.inputs = s['inputs']

        self.seeds = False
        self.evnts = False
        self.trange = False
        for i in self.inputs:
            if i['type'] == 'miniseed':
                self.miniseeds = [path.abspath(p['file']) for p in i['files']]
                log.info('Miniseeds were supplied for picking algorithm')
                self.seeds = True

            if i['type'] == 'events':
                self.events = i['events']
                log.info('Events were supplied for picking algorithm')
                self.evnts = True

            if i['type'] == 'time':
                self.time_range = i['times']
                self.time_range = {k: UTCDateTime(v) for k, v
                                   in self.time_range.items()}
                log.info('Time range was supplied for picking algorithm')
                self.trange = True

        sum_inputs = self.seeds + self.evnts + self.trange
        if sum_inputs != 1:
            raise ConfigException('Only one of miniseed, events or time '
                                  'range can be specified')

        self.picker = s['picker']

        self.detrend = s['detrend'] if 'detrend' in s else False

        self.filter = s['filter'] if 'filter' in s else False


class ConfigException(Exception):
    pass
