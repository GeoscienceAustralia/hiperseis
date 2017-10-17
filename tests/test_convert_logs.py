import sys
import os
from os.path import join, splitext
import json
from subprocess import check_call
import pytest

TESTS = os.path.dirname(__file__)
PASSIVE = os.path.dirname(TESTS)
LOGS = join(PASSIVE, 'convert_logs')
TESTDATA = join(LOGS, 'testdata')

test_dat_files = ['SQ2F2LOGFILE.dat', 'AQ3G4LOGFILE.dat']
PY3 = (sys.version_info[0] == 3)


@pytest.fixture(params=test_dat_files)
def logfile(request):
    return request.param


@pytest.mark.skipif(PY3, reason='only supported in py27')
def test_convert_logs(random_filename, logfile):
    output_dir = os.path.dirname(random_filename())

    check_call(('python', join(LOGS, 'decode_datfile.py'),
                    join(TESTDATA, logfile), '-o', output_dir,
                    '-y', '2015'))

    basename = splitext(logfile)[0]
    d_out = json.load(open(join(output_dir, basename + '.json'), 'r'))
    d_stored = json.load(open(join(TESTDATA, basename + '.json'), 'r'))
    assert d_out == d_stored


@pytest.mark.skipif(PY3, reason='only supported in py27')
def test_convert_logs_parallel(random_filename):
    output_dir = os.path.dirname(random_filename())
    check_call(('python', join(LOGS, 'decode_datfile.py'),
                TESTDATA, '-o', output_dir,
                '-y', '2015'))
    for d in test_dat_files:
        basename = splitext(d)[0]
        d_out = json.load(open(join(output_dir, basename + '.json'), 'r'))
        d_stored = json.load(open(join(TESTDATA, basename + '.json'), 'r'))
        print(open(join(TESTDATA, basename + '.json'), 'r'),
              open(join(output_dir, basename + '.json'), 'r'))
        assert d_out == d_stored
