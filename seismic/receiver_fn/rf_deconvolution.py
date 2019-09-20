#!/usr/bin/env python

import numpy as np
import scipy
import scipy.signal
import matplotlib.pyplot as plt

import obspy
import rf

from seismic.receiver_fn.rf_h5_file_event_iterator import IterRfH5FileEvents

# pylint: disable=invalid-name


def fcorrelate(f, g):
    """
    correlation routine - correlates f and g, normalized
      by the zero-lag autocorrelation of g. Returns only
      part of correlation corresponding to positive lag.
    """
    n = len(f)
    assert len(g) == n

    # Use the Numerical Recipes routine to compute the cross correlation
    #     call correl(f, g, n2, c) # output ends up in c
    c = scipy.signal.correlate(f, g)
    # Zero or positive lags in result are in elements [-n:] of c.
    c_pos = c[-n:]

    # compute the zero-lag autocorrelation of g
    sum0 = np.dot(g, g)
    normalize_factor = 1.0 / sum0
    xc = c_pos * normalize_factor

    return xc


def get_residual(x, y):
    """
    get the residual between x and y

    n = len(x), len(y)
    """
    r = x - y
    sumsq = np.dot(r, r)

    return r, sumsq


def gauss_filter(x, gwidth_factor, dt):

    """
    Convolve a function with a unit-area Gaussian filter.
    i.e. apply low-pass filter to input x using Gaussian function in freq domain.
    """
    n = len(x)
    two_pi = 2 * np.pi

    fft_x = np.fft.rfft(x)  # complex array
    n2 = len(fft_x)

    df = 1.0 / (float(n) * dt)
    d_omega = two_pi * df
    gwidth = 4.0 * gwidth_factor * gwidth_factor

    omega = np.arange(n2) * d_omega
    gauss = np.exp(-omega * omega / gwidth)
    fft_x = fft_x * gauss

    x_filt = np.fft.irfft(fft_x, n)  # real_array

    return x_filt


def phase_shift(x, time_shift, dt):
    """
    phase shifts a signal
    """
    n = len(x)
    two_pi = 2 * np.pi

    fft_x = np.fft.rfft(x)  # complex array
    n2 = len(fft_x)

    df = 1 / (float(n) * dt)
    d_omega = two_pi * df

    omega = np.arange(n2) * d_omega
    # The phase shift is omega (angular velocity) * delta_t (the time shift).
    spectral_shift = np.exp(omega*time_shift*1j)
    fft_x = fft_x * spectral_shift

    x_shifted = np.fft.irfft(fft_x, n)  # real_array

    return x_shifted


def build_decon(amps, shifts, n, gwidth, dt):
    """
    compute the predicted time series from a set of
    amplitudes and shifts
    """
    p = np.zeros(n)
    for i, amp in zip(shifts, amps):
        p[i] += amp

    p_filt = gauss_filter(p, gwidth, dt)

    return p_filt, p


def _convolve(x, y):
    """
    """
    n = len(x)
    assert len(y) == n
    z = scipy.signal.convolve(x, y)
    z_sync = z[:n]
    return z_sync


def iterdeconfd(denominator, numerator, maxbumps, tol=1.0e-3, gwidth=2.0, time_shift=None, lpositive=False):
    """
    Iterative deconvolution of source and response signal to generate seismic receiver function.

    :param denominator: The source signal
    :param response: The response signal (e.g. R or T component)
    :return:
    """
    MAXPTS = 100000
    assert len(denominator) <= MAXPTS  # Sanity check input size
    assert len(numerator) <= MAXPTS
    MAX_PULSES = 1000  # Maximum number of pulses to synthesize in p

    amps = []
    shifts = []

    print()
    print('Program iterdeconfd - Version 1.04, 1997-98')
    print('Chuck Ammon, Saint Louis University')
    print('Adapted to Python by Andrew Medlin, Geoscience Australia (2019)')
    print()

    if maxbumps > MAX_PULSES:
        print('Maximum Number of bumps is %d' % MAX_PULSES)
        maxbumps = MAX_PULSES
    # end if

    if time_shift is None:
        try:
            time_shift = float(denominator.stats.onset - denominator.stats.starttime)
        except AttributeError:
            time_shift = 0
    # end if

    # Clone numerator data
    f = numerator.copy().data
    npts = len(f)

    # Clone denominator data
    g = denominator.copy().data
    assert len(g) == npts, "g should be same length as f"

    dt = 1.0 / denominator.meta.sampling_rate

    if npts > MAXPTS:
        print('Too many points needed.')
        print('npts = %d' % npts)
        exit(1)
    # end if

    # Now begin the cross-correlation procedure

    # Put the filter in the signals
    f_hat = gauss_filter(f, gwidth, dt)
    g_hat = gauss_filter(g, gwidth, dt)

    # compute the power in the "numerator" for error scaling
    power = np.dot(f_hat, f_hat)

    # correlate the signals
    xc = fcorrelate(f_hat, g_hat)

    # find the peak in the correlation
    # maxlag = npts // 2
    # print('The maximum spike delay is %g sec' % (float(maxlag) * dt))

    if lpositive:
        shifts.append(np.argmax(xc))
    else:
        xc_abs = np.abs(xc)
        shifts.append(np.argmax(xc_abs))
    # end if
    amps.append(xc[shifts[-1]])

    nshifts = 0

    # compute the predicted deconvolution result
    p, pulses = build_decon(amps, shifts, npts, gwidth, dt)

    # convolve the prediction with the denominator signal
    f_hat_predicted = _convolve(p, g)

    # compute the residual (initial error is 1.0)
    resid, sumsq_ip1 = get_residual(f_hat, f_hat_predicted)

    sumsq_i = 1.0
    sumsq_ip1 = sumsq_ip1 / power
    d_error = 100 * (sumsq_i - sumsq_ip1)

    print('%11s %s' % ('Iteration', '  Spike amplitude  Spike delay   Misfit   Improvement'))
    print('%10d  %16.9e  %10.3f   %7.2f%%   %9.4f%%' % (
        nshifts, dt * amps[-1], (shifts[-1] - 1) * dt, 100 * sumsq_ip1, d_error))

    # *************************************************************************
    while (np.abs(d_error) > tol) and (nshifts < maxbumps):

        nshifts = nshifts + 1
        sumsq_i = sumsq_ip1

        xc = fcorrelate(resid, g_hat)

        if lpositive:
            shifts.append(np.argmax(xc))
        else:
            xc_abs = np.abs(xc)
            shifts.append(np.argmax(xc_abs))
        # end if
        amps.append(xc[shifts[-1]])

        p, pulses = build_decon(amps, shifts, npts, gwidth, dt)

        f_hat_predicted = _convolve(p, g)

        resid, sumsq_ip1 = get_residual(f_hat, f_hat_predicted)
        sumsq_ip1 = sumsq_ip1 / power
        d_error = 100 * (sumsq_i - sumsq_ip1)

        print('%10d  %16.9e  %10.3f   %7.2f%%   %9.4f%%' % (nshifts,
                                                            dt * amps[-1], shifts[-1] * dt,
                                                            100 * sumsq_ip1, d_error))
    # end while

    # *************************************************************************

    print('Last Error Change = %9.4f%%' % d_error)

    # if the last change made no difference, drop it
    fit = 100 - 100 * sumsq_ip1
    if d_error <= tol:
        nshifts = nshifts - 1
        fit = 100 - 100 * sumsq_i
        print('Hit the min improvement tolerance - halting.')
    # end if

    if nshifts >= maxbumps:
        print('Hit the max number of bumps - halting.')
    # end if

    print('Number of bumps in final result: %d' % nshifts)
    print('The final deconvolution reproduces %6.1f%% of the signal.' % fit)

    # *************************************************************************

    # compute the final prediction
    p, pulses = build_decon(amps, shifts, npts, gwidth, dt)

    response_trace = denominator.copy()  # clone input trace
    f_hat_predicted = _convolve(p, g)
    response_trace.data = f_hat_predicted

    p_shifted = phase_shift(p, -time_shift, dt)
    rf_trace = numerator.copy()  # clone input trace
    rf_trace.data = p_shifted

    return rf_trace, response_trace, pulses

# end func


def iter3c(stream):
    return rf.util.IterMultipleComponents(stream, key='onset', number_components=(2, 3))

#------------------------------------------------------------------------------

def main():
    """
    Main
    """
    use_Ammon_data = False
    if use_Ammon_data:
        src_sac_z = r"C:\software\hiperseis\seismic\receiver_fn\DATA\Ammon_test\lac_sp.z"
        src_sac_r = r"C:\software\hiperseis\seismic\receiver_fn\DATA\Ammon_test\lac_sp.r"
        src_expected_rf = r"C:\software\hiperseis\seismic\receiver_fn\DATA\Ammon_test\lac.i.eqr"
        source = obspy.read(src_sac_z, format='sac')
        response = obspy.read(src_sac_r, format='sac')
        expected_rf = obspy.read(src_expected_rf, format='sac')
        rf_trace, predicted_response, pulses = iterdeconfd(source[0], response[0], 200, gwidth=2.5, time_shift=5.0)
        expected_rf.plot()
        rf_trace.plot()
        predicted_response.plot()
        plt.stem(pulses)
        plt.show()
    else:
        # Load test data
        src_trace_file = r"C:\software\hiperseis\seismic\receiver_fn\DATA\7W.BL05_event_waveforms_for_rf_filtered.h5"
        src_data = obspy.read(src_trace_file, format='h5')
        # Run deconv on traces associated with same events
        max_iterations = 200
        time_window = (-10, 30)
        for stream3c in iter3c(src_data.copy()):
            stream3c.filter('lowpass', freq=3.0, corners=2, zerophase=True).interpolate(10.0)
            rf_stream = rf.RFStream(stream3c)
            rf_stream.rotate('NE->RT')
            rf_stream.trim2(*time_window, reftime='onset')
            rf_stream.detrend('linear')
            rf_stream.taper(0.2, max_length=0.5)
            source = rf_stream .select(component='Z')[0]
            response = rf_stream .select(component='R')[0]
            rf_trace, predicted_response, pulses = iterdeconfd(source, response, max_iterations, gwidth=2.5)
            rf_trace.plot()
        # end for
    # end if

# end main

if __name__ == "__main__":
    main()
