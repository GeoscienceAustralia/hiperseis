from __future__ import division

import random
import string
import os.path
import pytest

TESTS = os.path.dirname(__file__)
DATA = os.path.join(TESTS, 'mocks', 'data')


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
