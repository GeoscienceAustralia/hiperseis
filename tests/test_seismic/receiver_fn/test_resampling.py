#!/usr/bin/env python
"""Unit testing for seismic resampling functions
"""

import pytest

import numpy as np
import scipy as sp

from seismic.stream_processing import sinc_resampling


def test_sinc_resampling():
    # Test upsampling x2
    x = sp.signal.ricker(16, 1.5)
    t = np.linspace(0, 5, 16)
    t_new = np.linspace(0, 5, 2*16 - 1)
    x_new = sinc_resampling(t, x, t_new)
    assert np.allclose(x, x_new[0:32:2])

    # Test upsampling to arbitrary, non-uniform new time values.
    # Expect the interpolated values to be close to what we get
    # using linear interpolation.
    t = np.linspace(0, 5, 17)
    x = sp.signal.gausspulse(t - t.mean(), fc=4, bw=0.12)
    t_new = np.linspace(0, 5, 51)
    nan_idx = np.random.randint(1, len(t_new) - 2, 10)
    t_new[nan_idx] = np.NaN
    t_new = t_new[~np.isnan(t_new)]
    x_new = sinc_resampling(t, x, t_new)
    # Perform high frequency, high order interpolation and check values are close.
    t_hf = np.linspace(0, 5, 1000)
    x_hf = sp.interpolate.interp1d(t, x)(t_hf)
    x_new_hf = sp.interpolate.interp1d(t_new, x_new)(t_hf)
    assert np.sqrt(np.mean((x_new_hf - x_hf)**2)) < 0.1
# end func


if __name__ == "__main__":
    test_sinc_resampling()
# end if
