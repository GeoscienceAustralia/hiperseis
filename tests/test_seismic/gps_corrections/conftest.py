#!/usr/bin/env python

import os
import pytest

from seismic.gps_corrections import picks_reader_utils as pru


def _read_picks():
    test_ensemble_file = os.path.join(os.path.split(__file__)[0], 'data', 'test_picks_ensemble.txt')
    picks = pru.read_picks_ensemble(test_ensemble_file)
    return picks


def _iris_stntxt_file():
    test_iris_stntxt_file = os.path.join(os.path.split(__file__)[0], 'data', 'test_AU_GE_IRIS_inv.txt')
    return test_iris_stntxt_file


@pytest.fixture(scope='session')
def df_picks():
    return _read_picks()


@pytest.fixture(scope='session')
def iris_stntxt_file():
    return _iris_stntxt_file()
