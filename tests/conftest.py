from __future__ import division
import yaml
import random
import string
import os.path
import pytest
import datetime

from seismic import config

TESTS = os.path.dirname(__file__)
DATA = os.path.join(TESTS, 'mocks', 'data')

TESTS = os.path.dirname(__file__)
MOCKS = os.path.join(TESTS, 'mocks')


@pytest.fixture
def mseed():
    return os.path.join(MOCKS, 'ga2017qxlpiu_short.mseed')


@pytest.fixture
def xml():
    return os.path.join(MOCKS, 'events', 'ga2017qxlpiu.xml')


@pytest.fixture
def test_dir():
    return TESTS


@pytest.fixture
def data_dir():
    return DATA


@pytest.fixture
def random_filename(tmpdir_factory):
    def make_random_filename(ext=''):
        dir = str(tmpdir_factory.mktemp('seismic').realpath())
        fname = ''.join(random.choice(string.ascii_lowercase)
                        for _ in range(10))
        return os.path.join(dir, fname + ext)
    return make_random_filename


@pytest.fixture
def conf():
    return {'inputs': [],
            'picker': None}


@pytest.fixture
def events():
    return {'name': 'my events',
            'events': [{'event_id': 'ga2017abcdefg'},
                       {'event_id': 'ga2017hijklmn'}],
            'type': 'events'}


@pytest.fixture
def times():
    return {'times': {'end_time': datetime.datetime(2017, 3, 28, 17, 18, 30),
                      'start_time':
                          datetime.datetime(2017, 3, 28, 16, 18, 30)},
            'name': 'my time range', 'type': 'time'}


@pytest.fixture
def miniseed():
    return {'files': [{'file': os.path.join(DATA, 'ev0_6.a01.gse2')},
                      {'file': os.path.join(DATA, 'ev0_6.a02.gse2')}],
            'name': 'my miniseed files', 'type': 'miniseed'}


@pytest.fixture
def params_dict(events, times, miniseed):
    return {
        'miniseeds': miniseed,
        'events': events,
        'times': times
        }


@pytest.fixture
def bandpass_filter():
    return {'type': 'bandpass',
            'params': {'freqmin': 2.0,
                       'freqmax': 16.0,
                       'corners': 3,
                       'zerophase': True}}


@pytest.fixture
def miniseed_conf(random_filename, conf, params_dict, bandpass_filter):
    conf['inputs'].append(params_dict['miniseeds'])
    conf['filter'] = bandpass_filter
    yaml_file = random_filename(ext='.yml')
    with open(yaml_file, 'w') as outfile:
        yaml.dump(conf, outfile, default_flow_style=False)
    cf = config.Config(yaml_file)
    return cf
