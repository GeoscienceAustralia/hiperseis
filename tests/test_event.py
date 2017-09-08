import os
from obspy import read_events


def test_obspy_read_sc3ml(test_dir):
    read_events(os.path.join(test_dir, 'mocks', 'ga2017qxlpiu.xml'),
                format='SC3ML')
