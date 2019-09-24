#!/usr/bin/env python

import numpy as np
from scipy import signal

# import matplotlib.pyplot as plt

from seismic.receiver_fn.rf_synthetic import synthesize_rf_dataset
from seismic.receiver_fn.rf_deconvolution import iter_deconv_pulsetrain

# pylint: disable=invalid-name, missing-docstring


def test_ammon_iter_deconv():
    # Synthesize known RF as a series of delta functions and verify that iterative
    # deconvolution recovers original RF exactly, to within numerical tolerance.
    inclinations = np.array([20.0])
    distances = np.array([60.0])
    amplitudes = [1, 0.4, 0.2]
    H = 42.0  # km
    V_p = 6.4  # km/s
    V_s = 3.8  # km/s
    F_s = 10.0  # Hz
    # This function call generates the container for the synthetic RF, then we overwrite it with pure delta functions
    rf_radial, arrivals = synthesize_rf_dataset(H, V_p, V_s, inclinations, distances, F_s, amplitudes=amplitudes)
    rf_radial = rf_radial[0]
    # Patch over the filtered RF data with pure delta functions for testing purposes
    times = rf_radial.times() - (rf_radial.stats.onset - rf_radial.stats.starttime)
    rf_radial.data = np.zeros_like(rf_radial.data)
    rf_radial.data[np.round((np.array(arrivals) - times[0])*F_s).astype(int)] = amplitudes
    rf_radial.data /= F_s
    # Generate synthetic vertical seismic trace as a filtered decaynig cosine. The exact form here shouldn't matter,
    # but better that it be plausible.
    g = 1000*np.exp(-times/4)*np.cos(2*np.pi*1*times)
    g[times < 0] = 0
    # Apply weak filter to g to smooth onset, cutting off above 2.0 Hz
    bw_filt = signal.butter(2, 2.0/F_s)
    g = signal.filtfilt(bw_filt[0], bw_filt[1], g)

    ## Test scaling. Deconvolving a function with itself should return pulse of unit area with maximum at onset.
    denominator = rf_radial.copy()
    denominator.data = g
    numerator = rf_radial.copy()
    numerator.data = g
    max_pulses = 50
    rf_trace, pulses, f_hat, response_trace, fit = iter_deconv_pulsetrain(denominator, numerator, max_pulses)
    assert np.isclose(np.sum(rf_trace), 1.0)
    assert times[np.nonzero(rf_trace.data == np.max(rf_trace.data))][0] == 0.0

    # Generate radial seismic trace corresponding to source signal g convolved with the RF.
    f = signal.convolve(rf_radial.data, g, mode='full')
    n_leadin = len(times[times < 0])
    f = f[n_leadin:n_leadin+len(times)]

    # Apply iterative deconvolution
    denominator.data = g
    numerator.data = f
    rf_trace, pulses, f_hat, response_trace, fit = iter_deconv_pulsetrain(denominator, numerator, max_pulses)
    # Check that predicted response is very close to filtered observation
    assert np.allclose(f_hat, response_trace, rtol=2.0e-3, atol=1.0e-3)
    # Check that we have the exact number of delta function pulses that we synthesized
    assert np.sum(pulses > 0) == 3
    # Check that the fit is perfect, which it should be when RF consists of perfect delta functions.
    assert np.isclose(100.0, fit)
    # Check that amplitude of the computed RF matches the synthesized RF.
    assert np.isclose(np.sum(amplitudes)/F_s, np.sum(rf_trace), rtol=1.0e-3)

# end func


if __name__ == "__main__":
    test_ammon_iter_deconv()
