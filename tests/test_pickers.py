import os
import pytest
from obspy.core import read as obspy_read
from seismic.pickers import pickermaps

TESTS = os.path.dirname(__file__)
MOCKS = os.path.join(TESTS, 'mocks')
mseed = os.path.join(MOCKS, 'ga2017qxlpiu_short.mseed')
xml = os.path.join(MOCKS, 'ga2017qxlpiu.xml')

algos = list(pickermaps.keys())


@pytest.fixture(params=algos)
def algorithm(request):
    return request.param


def test_pickermaps(algorithm):
    """
    basic operation test for now
    """
    picker = pickermaps[algorithm]()
    st = obspy_read(mseed)
    for s in st[:2]:
        picker.picks(s)
