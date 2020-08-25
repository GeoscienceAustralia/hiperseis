import os

import pytest
import numpy as np

from seismic.receiver_fn.ccp_correction import correct
from tests.conftest import TESTS


@pytest.fixture()
def data_dir():
    return os.path.join(TESTS, 'test_seismic', 'receiver_fn', 'data')


def test_ccp_correction(tmpdir, data_dir):
    ccp_data = os.path.join(data_dir, 'ccp_test.csv')
    corr_data = os.path.join(data_dir, 'hk_test.csv')
    outfile = os.path.join(tmpdir, 'test.csv')
    correct(ccp_data, corr_data, outfile)
    test_result = np.genfromtxt(outfile, delimiter=',', dtype=None, encoding=None,
                                names=['sta', 'lon', 'lat', 'depth', 'weight'])
    assert test_result.shape == (7,)
    expected_depths = np.array([2, 7, 12, 30, 30, 35, 50])
    np.testing.assert_equal(test_result['depth'], expected_depths)
