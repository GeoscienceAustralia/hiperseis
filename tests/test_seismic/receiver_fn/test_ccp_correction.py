import os

import pytest
import numpy as np

from seismic.receiver_fn.ccp_correction import correct


def test_ccp_correction(tmpdir):
    ccp_data = 'data/ccp_test.csv'
    corr_data= 'data/hk_test.csv'
    outfile = os.path.join(tmpdir, 'test.csv')
    correct(ccp_data, corr_data, outfile)
    test_result = np.genfromtxt(outfile, delimiter=',', dtype=None, encoding=None,
                                names=['sta', 'lon', 'lat', 'depth', 'weight'])
    assert test_result.shape == (7,)
    expected_depths = np.array([2, 7, 12, 30, 30, 35, 50])
    np.testing.assert_equal(test_result['depth'], expected_depths)
