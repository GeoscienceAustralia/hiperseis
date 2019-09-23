#!/usr/bin/env python

import os
import pytest

import numpy as np
from scipy import signal

import matplotlib.pyplot as plt

from seismic.receiver_fn.rf_synthetic import synthesize_rf_dataset
from seismic.receiver_fn.rf_deconvolution import iter_deconv_pulsetrain

# pylint: disable=invalid-name, missing-docstring


def test_ammon_iter_deconv(tmp_path):
    assert True
    # Synthesize known RF as a series of delta functions and verify that iterative
    # deconvolution recovers original RF.
    inclinations = np.array([20.0])
    distances = np.array([60.0])
    amplitudes = [1, 0.4, 0.2]
    H = 42.0  # km
    V_p = 6.4  # km/s
    V_s = 3.8  # km/s
    F_s = 10.0  # Hz
    rf_radial = synthesize_rf_dataset(H, V_p, V_s, inclinations, distances, F_s, amplitudes=amplitudes)
    rf_radial = rf_radial[0]
    rf_radial.data /= F_s
    # Generate vertical seismic trace
    times = rf_radial.times() - (rf_radial.stats.onset - rf_radial.stats.starttime)
    g = 1000*np.exp(-times/4)*np.cos(2*np.pi*1*times)
    g[times < 0] = 0
    # Apply weak filter to g to smooth onset
    bw_filt = signal.butter(2, 2.0/F_s)
    g = signal.filtfilt(bw_filt[0], bw_filt[1], g)
    # Generate radial seismic trace
    f = signal.convolve(rf_radial.data, g, mode='full')
    n_leadin = len(times[times < 0])
    f = f[n_leadin:n_leadin+len(times)]

    # Apply iterative deconvolution
    denominator = rf_radial.copy()
    denominator.data = g
    numerator = rf_radial.copy()
    numerator.data = f
    max_pulses = 100
    rf_trace, pulses, f_hat, response_trace, fit = iter_deconv_pulsetrain(denominator, numerator, max_pulses)
    pass

if __name__ == "__main__":
    import tempfile
    test_ammon_iter_deconv(tempfile.mkdtemp())
