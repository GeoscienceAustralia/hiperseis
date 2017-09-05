import logging
import yaml
from os import path

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
            s = yaml.load(f)

        self.name = path.basename(yaml_file).rsplit(".", 1)[0]
        self.inputs = s['inputs']

        for i in self.inputs:
            if i['type'] == 'miniseed':
                self.miniseeds = [path.abspath(p['file']) for p in i['files']]

            if i['type'] == 'events':
                self.events = i['events']
