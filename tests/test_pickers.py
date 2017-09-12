import os
import yaml
import pytest
from obspy.core import read as obspy_read
from seismic.pickers import pickermaps
from seismic import config

TESTS = os.path.dirname(__file__)
MOCKS = os.path.join(TESTS, 'mocks')
mseed = os.path.join(MOCKS, 'ga2017qxlpiu.mseed')
xml = os.path.join(MOCKS, 'ga2017qxlpiu.xml')

algos = list(pickermaps.keys())


@pytest.fixture(params=algos)
def algorithm(request):
    return request.param


def test_pickermaps(conf, miniseed, random_filename, algorithm):
    """
    basic operation test for now
    """
    conf['inputs'].append(miniseed)
    yaml_file = random_filename(ext='.yaml')
    with open(yaml_file, 'w') as outfile:
        yaml.dump(conf, outfile, default_flow_style=False)
    cf = config.Config(yaml_file)
    picker = pickermaps[algorithm]()
    st = obspy_read()
    for s in st[:2]:
        picker.picks(s)
