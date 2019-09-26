#!/usr/bin/env python
"""Custom receiver function deconvolution methods and support functions.
"""

import logging

import numpy as np
import scipy
import scipy.signal

# pylint: disable=invalid-name,too-many-locals

logging.basicConfig()


def _xcorrelate(f, g):
    """
    Correlation routine - correlates f and g, normalized by the zero-lag autocorrelation of g.
    Returns only part of correlation corresponding to positive lag.
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
# end func


def _get_residual(x, y):
    """
    Get the residual between x and y and the power of the residual.
    """
    r = x - y
    sumsq = np.dot(r, r)

    return r, sumsq
# end func


def _gauss_filter(x, gwidth_factor, dt):
    """
    Apply low-pass filter to input x using Gaussian function in freq domain.
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
# end func


def _build_decon(amps, shifts, n, gwidth, dt):
    """
    Compute the predicted RF time series from a set of pulse amplitudes and time shifts.
    """
    p = np.zeros(n)
    for i, amp in zip(shifts, amps):
        p[i] += amp

    p_filt = _gauss_filter(p, gwidth, dt)

    return p_filt, p
# end func


def _convolve(x, y, leadin):
    """
    Helper function for convolving approximate receiver function with source signal
    to get estimated observation.
    """
    n = len(x)
    assert len(y) == n
    z = scipy.signal.convolve(x, y, mode='full')
    z_sync = z[leadin:leadin + n]
    return z_sync
# end func


def iter_deconv_pulsetrain(numerator, denominator, sampling_rate, time_shift, max_pulses=1000,
                           tol=1.0e-3, gwidth=2.5, only_positive=False, log=None, **kwargs):
    """
    Iterative deconvolution of source and response signal to generate seismic receiver function.
    Adapted to Python by Andrew Medlin, Geoscience Australia (2019), from Chuck Ammon's
    (Saint Louis University) `iterdeconfd` Fortran code, version 1.04.

    Note this is not really a frequency-domain deconvolution method, since the deconvolution is
    based on generating a time-domain pulse train which is filtered and convolved with source
    signal in time domain in order to try to replicate observation. Results should be the same
    even if some of the spectral techniques used in functions (such as `gauss_filter`) were replaced
    by non-spectral equivalents.

    :param numerator: The observed response signal (e.g. R, Q or T component)
    :type numerator: numpy.array(float)
    :param denominator: The source signal (e.g. L or Z component)
    :type denominator: numpy.array(float)
    :param sampling_rate: The sampling rate in Hz of the numerator and denominator signals
    :type sampling_rate: float
    :param time_shift: Time shift (sec) from start of input signals until expected P-wave arrival onset.
    :type time_shift: float
    :param max_pulses: Maximum number of delta function pulses to synthesize in the unfiltered RF (up to 1000)
    :type max_pulses: int
    :param tol: Convergence tolerance, iteration stops if change in error falls below this value
    :type tol: float
    :param gwidth: Gaussian filter width in the normalized frequency domain.
    :type gwidth: float
    :param only_positive: If true, only use positive pulses in the RF.
    :type only_positive: bool
    :param log: Log instance to log messages to.
    :type log: logging.Logger
    :return: RF trace, pulse train, expected response signal, predicted response signal, quality of fit statistic
    :rtype: numpy.array(float), numpy.array(float), numpy.array(float), numpy.array(float), float
    """
    if log is None:
        log = logging.getLogger(__name__)
        log.setLevel(logging.INFO)
    # end if

    MAXPTS = 100000
    assert len(denominator) <= MAXPTS  # Sanity check input size
    assert len(numerator) <= MAXPTS
    MAX_PULSES = 1000  # Maximum number of pulses to synthesize in p

    amps = []
    shifts = []

    if max_pulses > MAX_PULSES:
        log.warning('Maximum Number of pulses is %d, clipping input value', MAX_PULSES)
        max_pulses = MAX_PULSES
    # end if

    # Clone numerator data
    f = numerator.copy()
    npts = len(f)

    # Clone denominator data
    g = denominator.copy()
    assert len(g) == npts, "g should be same length as f"

    dt = 1.0 / sampling_rate

    # Now begin the cross-correlation procedure

    # Put the filter in the signals
    f_hat = _gauss_filter(f, gwidth, dt)
    g_hat = _gauss_filter(g, gwidth, dt)

    # compute the power in the "numerator" for error scaling
    power = np.dot(f_hat, f_hat)

    # correlate the signals
    xc = _xcorrelate(f_hat, g_hat)

    # Exclude non-causal offsets from being added to the pulse locations.
    non_causal_leadin = int(np.ceil(time_shift/dt))
    max_causal_idx = npts - non_causal_leadin
    assert max_causal_idx > 0
    if only_positive:
        shift_index = np.argmax(xc[:max_causal_idx])
    else:
        xc_abs = np.abs(xc)
        shift_index = np.argmax(xc_abs[:max_causal_idx])
    # end if
    shifts.append(shift_index + non_causal_leadin)
    amps.append(xc[shift_index])

    num_pulses = 1

    # compute the predicted deconvolution result
    p, pulses = _build_decon(amps, shifts, npts, gwidth, dt)

    # convolve the prediction with the denominator signal
    f_hat_predicted = _convolve(p, g, non_causal_leadin)

    # compute the residual (initial error is 1.0)
    resid, sumsq_ip1 = _get_residual(f_hat, f_hat_predicted)

    sumsq_i = 1.0
    sumsq_ip1 = sumsq_ip1 / power
    d_error = 100 * (sumsq_i - sumsq_ip1)

    log.info('%11s %s', 'Iteration', '  Spike amplitude  Spike delay   Misfit   Improvement')
    log.info('%10d  %16.9e  %10.3f   %7.2f%%   %9.4f%%', num_pulses, dt * amps[-1], (shifts[-1] - 1) * dt,
             100 * sumsq_ip1, d_error)

    # *************************************************************************
    while (np.abs(d_error) > tol) and (num_pulses < max_pulses):

        num_pulses += 1
        sumsq_i = sumsq_ip1

        xc = _xcorrelate(resid, g_hat)

        if only_positive:
            shift_index = np.argmax(xc[:max_causal_idx])
        else:
            xc_abs = np.abs(xc)
            shift_index = np.argmax(xc_abs[:max_causal_idx])
        # end if
        shifts.append(shift_index + non_causal_leadin)
        amps.append(xc[shift_index])

        p, pulses = _build_decon(amps, shifts, npts, gwidth, dt)

        f_hat_predicted = _convolve(p, g, non_causal_leadin)

        resid, sumsq_ip1 = _get_residual(f_hat, f_hat_predicted)
        sumsq_ip1 = sumsq_ip1 / power
        d_error = 100 * (sumsq_i - sumsq_ip1)

        log.info('%10d  %16.9e  %10.3f   %7.2f%%   %9.4f%%', num_pulses, dt * amps[-1], shifts[-1] * dt,
                 100 * sumsq_ip1, d_error)
    # end while

    # *************************************************************************

    log.info('Last Error Change = %9.4f%%', d_error)

    # Compute final fit
    fit = 100 - 100 * sumsq_ip1

    if (num_pulses >= max_pulses) and (np.abs(d_error) > tol):
        log.warning('Hit the max number of pulses - not halting due to convergence.')
    # end if

    num_unique_pulses = len(set(shifts))
    log.info('Number of pulses in final result: %d', num_unique_pulses)
    log.info('The final deconvolution reproduces %5.1f%% of the signal.', fit)

    # *************************************************************************

    # compute the final prediction
    p, pulses = _build_decon(amps, shifts, npts, gwidth, dt)

    f_hat_predicted = _convolve(p, g, non_causal_leadin)
    rf_trace = p

    return rf_trace, pulses, f_hat, f_hat_predicted, fit

# end func


def rf_iter_deconv(response_data, source_data, sr, tshift, min_fit_warning=80.0, normalize=False, **kwargs):
    """Adapter function to rf library.  To use, add arguments `deconvolve='func', func=rf_iter_deconv` to
    `rf.RFStream.rf()` function call.

    :param response_data: List of response signals for which to compute receiver function
    :type response_data: list of numpy.array(float)
    :param source_data: Source signal to use for computing the receiver functions
    :type source_data: numpy.array(float)
    :param sampling_rate: Sampling rate of the input signals (Hz)
    :type sampling_rate: float
    :param time_shift: Time shift (seconds) from start of signal to onset
    :type time_shift: float
    :param min_fit_warning: Minimum percentage of fit to include trace in results,
        otherwise will be returned as empty numpy array.
    :type min_fit_warning: float
    :param normalize: Flag indicating whether RF should be normalized. If True, will be scaled
        so that sqrt of sum of squares is 1. If RF is normalized, it will not be able to be used
        to directly generate response signal from source signal.
    :return: Receiver functions corresponding to the list of input response signals.
    :rtype: list of numpy.array(float)
    """
    normalize = bool(normalize)
    sampling_rate = sr
    time_shift = tshift
    denominator = source_data
    receiver_fns = []
    log = logging.getLogger(__name__)
    for numerator in response_data:
        # Any non-default parameters of deconvolution should be packed into kwargs
        rf_trace, _, _, _, fit = iter_deconv_pulsetrain(numerator, denominator, sampling_rate, time_shift, **kwargs)
        if fit < min_fit_warning:
            receiver_fns.append(np.array())
            log.warning("RF fit {:.2f}% below minimum acceptable threshold, discarding".format(fit))
        else:
            if normalize:
                sum_sq = np.sum(np.square(rf_trace))
                rf_trace /= np.sqrt(sum_sq)
            # end if
            receiver_fns.append(rf_trace)
        # end if
    # end for

    return receiver_fns
# end func
