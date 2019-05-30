#!/usr/bin/env python

import os
import pytest

from seismic.gps_corrections import picks_reader_utils as pru

def _read_picks():
    test_ensemble_file = os.path.join(os.path.split(__file__)[0], 'data', 'test_picks_ensemble.txt')
    return pru.read_picks_ensemble(test_ensemble_file)

@pytest.fixture(scope='session')
def df_picks():
    return _read_picks()
