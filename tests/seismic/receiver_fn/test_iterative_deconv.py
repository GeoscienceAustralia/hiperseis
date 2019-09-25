#!/usr/bin/env python
"""Unit testing for iterative deconvolution
"""

import numpy as np
from scipy import signal

from seismic.receiver_fn.rf_synthetic import synthesize_rf_dataset
from seismic.receiver_fn.rf_deconvolution import iter_deconv_pulsetrain, rf_iter_deconv

# pylint: disable=invalid-name, missing-docstring


def _generate_rf_radial(inclinations, distances, amplitudes):
    # Synthesize known RF as a series of delta functions
    assert len(distances) == len(inclinations)
    H = 42.0  # km
    V_p = 6.4  # km/s
    V_s = 3.8  # km/s
    F_s = 10.0  # Hz
    # This function call generates the container for the synthetic RF, then we overwrite it with pure delta functions
    rf_radial, arrivals = synthesize_rf_dataset(H, V_p, V_s, inclinations, distances, F_s, amplitudes=amplitudes)
    for tr in rf_radial:
        # Patch over the filtered RF data with pure delta functions for testing purposes
        time_shift = tr.stats.onset - tr.stats.starttime
        times = tr.times() - time_shift
        tr.data[:] = 0
        tr.data[np.round((np.array(arrivals) - times[0])*F_s).astype(int)] = amplitudes
        tr.data /= F_s
    # end for

    return rf_radial, F_s
# end func


def _generate_synthetic_source(times, amplitude=1000.0):
    # Generate synthetic vertical seismic trace as a filtered decaying cosine. The exact form here shouldn't matter,
    # but better that it be plausible.
    g = amplitude*np.exp(-times/4)*np.cos(2*np.pi*1*times)
    g[times < 0] = 0
    # Apply weak filter to g to smooth onset, cutting off above 2.0 Hz
    filter_order = 2
    cutoff = 0.2  # normalized cutoff freq
    bw_filt = signal.butter(filter_order, cutoff)
    g = signal.filtfilt(bw_filt[0], bw_filt[1], g)

    return g
# end func


def _generate_radial_from_src(pure_rf, source_signal, times):
    # Use Rf and source signal to generate radial test signal
    f = signal.convolve(pure_rf, source_signal, mode='full')
    n_leadin = len(times[times < 0])
    f = f[n_leadin:n_leadin+len(times)]

    return f
# end func


def test_ammon_iter_deconv():
    # Synthesize known RF and verify that iterative deconvolution recovers original RF exactly,
    # to within numerical tolerance.
    inclinations = np.array([20.0])
    distances = np.array([60.0])
    amplitudes = [1, 0.4, 0.2]
    rf_radial, F_s = _generate_rf_radial(inclinations, distances, amplitudes)
    rf_radial = rf_radial[0]
    time_shift = rf_radial.stats.onset - rf_radial.stats.starttime
    times = rf_radial.times() - time_shift

    # Generate synthetic vertical seismic trace.
    g = _generate_synthetic_source(times)

    ## Test scaling. Deconvolving a function with itself should return pulse of unit area with maximum at onset.
    numerator = rf_radial.copy().data
    denominator = rf_radial.copy().data
    assert F_s == rf_radial.stats.sampling_rate
    sampling_rate = F_s
    max_pulses = 50
    rf_trace, pulses, f_hat, response_trace, fit = iter_deconv_pulsetrain(numerator, denominator, sampling_rate,
                                                                          time_shift, max_pulses)
    assert np.isclose(np.sum(rf_trace), 1.0)
    assert times[np.nonzero(rf_trace == np.max(rf_trace))][0] == 0.0

    # Generate radial seismic trace corresponding to source signal g convolved with the RF.
    f = _generate_radial_from_src(rf_radial.data, g, times)

    # Apply iterative deconvolution
    numerator = f
    denominator = g
    rf_trace, pulses, f_hat, response_trace, fit = iter_deconv_pulsetrain(numerator, denominator, sampling_rate,
                                                                          time_shift, max_pulses)
    # Check that predicted response is very close to filtered observation
    assert np.allclose(f_hat, response_trace, rtol=2.0e-3, atol=1.0e-3)
    # Check that we have the exact number of delta function pulses that we synthesized
    assert np.sum(pulses > 0) == 3
    # Check that the fit is perfect, which it should be when RF consists of perfect delta functions.
    assert np.isclose(100.0, fit)
    # Check that amplitude of the computed RF matches the synthesized RF.
    assert np.isclose(np.sum(amplitudes)/F_s, np.sum(rf_trace), rtol=1.0e-3)

# end func


def test_rf_integration():
    import rf

    # Synthesize known RF.
    inclinations = np.array([20.0, 15.0, 10.0])
    distances = np.array([60.0, 70.0, 80.0])
    amplitudes = [1, 0.4, 0.2]
    rf_radial, F_s = _generate_rf_radial(inclinations, distances, amplitudes)

    time_shift = rf_radial[0].stats.onset - rf_radial[0].stats.starttime
    times = rf_radial[0].times() - time_shift

    # Generate synthetic vertical seismic trace.
    g = _generate_synthetic_source(times)

    # Collect test stream data with R and Z components into a RFStream object.
    rf_stream = rf.RFStream()
    np.random.seed(20190925)
    g_noise_scale = 5.0e-3 * np.abs(g).max()
    g_funcs = []
    f_funcs = []
    for i, tr in enumerate(rf_radial):
        # Assign trackable event id
        tr.stats.event_id = i
        # Add some random noise to the source signal
        g_noisy = g + np.random.normal(scale=g_noise_scale, size=g.shape)
        g_funcs.append(g_noisy)
        # Create source RFTrace
        src_tr = tr.copy()
        src_tr.data = g_noisy.copy()
        src_tr.stats.channel = 'HHZ'
        # Synthesize response signal and add some noise
        f = _generate_radial_from_src(tr.data, g_noisy, times)
        f_noise_scale = 5.0e-3 * np.abs(f).max()
        f_noisy = f + np.random.normal(scale=f_noise_scale, size=f.shape)
        f_funcs.append(f_noisy)
        # Create response RFTrace
        rsp_tr = tr.copy()
        rsp_tr.data = f_noisy.copy()
        rf_stream += src_tr
        rf_stream += rsp_tr
    # end for

    # Use rf library to compute comparative signals using time and freq domain deconvolution
    rf_freq = rf_stream.copy().rf(method='P', rotate=None, deconvolve='freq',
                                  response_components='R').select(component='R')
    try:
        rf_time = rf_stream.copy().rf(method='P', rotate=None, deconvolve='time',
                                      response_components='R').select(component='R')
    except NameError:
        import warnings
        # If Toeplitz not present on platform, rf may be unable to perform time-domain deconvolution
        warnings.warn("Unable to test default time-domain deconvolution from rf library")
        rf_time = None
    # end try

    # Call rf generator on rf.RFStream using our custom deconvolution function
    rf_iter = rf_stream.copy().rf(method='P', rotate=None, deconvolve='func', func=rf_iter_deconv,
                                  response_components='R', normalize=True).select(component='R')

    # Perform deconv directly and compare with rf_iter to check that rf call used our custom function.
    for i, (f, g) in enumerate(zip(f_funcs, g_funcs)):
        x, _, _, _, fit = iter_deconv_pulsetrain(f, g, F_s, time_shift)
        assert np.isclose(100.0, fit, rtol=1e-2)
        sum_sq = np.sum(np.square(x))
        x /= np.sqrt(sum_sq)
        assert rf_iter[i].stats.event_id == i
        assert np.allclose(rf_iter[i].data, x, rtol=1e-3, atol=5e-3)
        assert np.isclose(1.0, np.square(rf_iter[i].data).sum())
    # end for

    # Check that the local maxima of RF peaks found for different techniques all agree.
    # We expect exact agreement since the test data is very simple.
    def _local_maxima_mask_1d(arr):
        return (arr[1:-1] > arr[0:-2]) & (arr[1:-1] > arr[2:])
    # end func
    def _rms(arr):
        return np.sqrt(np.mean(np.square(arr)))
    # end func
    for i, tr in enumerate(rf_iter):
        d = tr.data
        # 1.3 factor here had to be tuned to exclude large amplitude wiggles in frequency domain deconvolution case
        expected_mask = _local_maxima_mask_1d(d) & (d[1:-1] > 1.3*_rms(d))
        mask_idx = np.nonzero(expected_mask)[0]
        d_fd = rf_freq[i].data
        fd_mask = _local_maxima_mask_1d(d_fd) & (d_fd[1:-1] > 1.3*_rms(d_fd))
        fd_mask_idx = np.nonzero(fd_mask)[0]
        assert np.all(fd_mask_idx == mask_idx)
        if rf_time is not None:
            d_td = rf_time[i].data
            td_mask = _local_maxima_mask_1d(d_td) & (d_td[1:-1] > 1.3*_rms(d_td))
            td_mask_idx = np.nonzero(td_mask)[0]
            assert np.all(td_mask_idx == mask_idx)
        # end if
    # end for

# end func

if __name__ == "__main__":
    test_ammon_iter_deconv()
    test_rf_integration()
