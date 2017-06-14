import os
from os.path import join
import json
from subprocess import check_call
import pytest

TESTS = os.path.dirname(__file__)
PASSIVE = os.path.dirname(TESTS)
LOGS = join(PASSIVE, 'convert_logs')
TESTDATA = join(LOGS, 'testdata')


@pytest.fixture(params=['SQ2F2LOGFILE', 'AQ3G4LOGFILE'])
def logfile(request):
    return request.param


def test_convert_logs(random_filename, logfile):
    output = random_filename(ext='.json')
    check_call(('python', join(LOGS, 'decode_datfile.py'),
                join(TESTDATA, logfile + '.dat'), '-o', output))
    d_out = json.load(open(output, 'r'))
    d_stored = json.load(open(join(TESTDATA, logfile + '.json'),
                              'r'))
    assert d_out == d_stored
