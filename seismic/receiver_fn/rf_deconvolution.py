#!/usr/bin/env python

import os
import logging

import numpy as np
import scipy
import scipy.signal
import matplotlib.pyplot as plt

import obspy
import rf

from seismic.receiver_fn.rf_plot_utils import plot_rf_stack

# pylint: disable=invalid-name

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


def get_residual(x, y):
    """
    Get the residual between x and y and the power of the residual.
    """
    r = x - y
    sumsq = np.dot(r, r)

    return r, sumsq


def gauss_filter(x, gwidth_factor, dt):
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


def phase_shift(x, time_shift, dt):
    """
    Phase shifts a signal using frequency domain technique for smooth arbitrary shifting.
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
    Compute the predicted RF time series from a set of pulse amplitudes and time shifts.
    """
    p = np.zeros(n)
    for i, amp in zip(shifts, amps):
        p[i] += amp

    p_filt = gauss_filter(p, gwidth, dt)

    return p_filt, p


def _convolve(x, y):
    """
    Helper function for convolving approximate receiver function with source signal
    to get estimated observation.
    """
    n = len(x)
    assert len(y) == n
    z = scipy.signal.convolve(x, y)
    z_sync = z[:n]
    return z_sync


def iter_deconv_pulsetrain(denominator, numerator, max_pulses, tol=1.0e-3, gwidth=2.5,
                           time_shift=None, lpositive=False):
    """
    Iterative deconvolution of source and response signal to generate seismic receiver function.
    Adapted to Python by Andrew Medlin, Geoscience Australia (2019), from Chuck Ammon's
    (Saint Louis University) `iterdeconfd` Fortran code, version 1.04.

    Note this is not really a frequency-domain deconvolution method, since the deconvolution is
    based on generating a time-domain pulse train which is filtered and convolved with source
    signal in time domain in order to try to replicate observation. Results should be the same
    even if some of the spectral techniques used in functions (such as `gauss_filter`) were replaced
    by non-spectral equivalents.

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

    if max_pulses > MAX_PULSES:
        log.warning('Maximum Number of pulses is %d, clipping input value', MAX_PULSES)
        max_pulses = MAX_PULSES
    # end if

    if time_shift is None:
        try:
            time_shift = float(denominator.stats.onset - denominator.stats.starttime)
        except AttributeError:
            time_shift = 0
    # end if
    assert time_shift is not None

    # Clone numerator data
    f = numerator.copy().data
    npts = len(f)

    # Clone denominator data
    g = denominator.copy().data
    assert len(g) == npts, "g should be same length as f"

    dt = 1.0 / denominator.meta.sampling_rate

    if npts > MAXPTS:
        log.error('Too many points needed.')
        log.error('npts = %d', npts)
        exit(1)
    # end if

    # Now begin the cross-correlation procedure

    # Put the filter in the signals
    f_hat = gauss_filter(f, gwidth, dt)
    g_hat = gauss_filter(g, gwidth, dt)

    # compute the power in the "numerator" for error scaling
    power = np.dot(f_hat, f_hat)

    # correlate the signals
    xc = _xcorrelate(f_hat, g_hat)

    # Exclude non-causal offsets from being added to the pulse locations.
    non_causal_leadin = int(np.ceil(time_shift/dt))
    max_causal_idx = npts - non_causal_leadin
    assert max_causal_idx > 0
    if lpositive:
        shifts.append(np.argmax(xc[:max_causal_idx]))
    else:
        xc_abs = np.abs(xc)
        shifts.append(np.argmax(xc_abs[:max_causal_idx]))
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

    log.info('%11s %s', 'Iteration', '  Spike amplitude  Spike delay   Misfit   Improvement')
    log.info('%10d  %16.9e  %10.3f   %7.2f%%   %9.4f%%', nshifts, dt * amps[-1], (shifts[-1] - 1) * dt,
             100 * sumsq_ip1, d_error)

    # *************************************************************************
    while (np.abs(d_error) > tol) and (nshifts < max_pulses):

        nshifts = nshifts + 1
        sumsq_i = sumsq_ip1

        xc = _xcorrelate(resid, g_hat)

        if lpositive:
            shifts.append(np.argmax(xc[:max_causal_idx]))
        else:
            xc_abs = np.abs(xc)
            shifts.append(np.argmax(xc_abs[:max_causal_idx]))
        # end if
        amps.append(xc[shifts[-1]])

        p, pulses = build_decon(amps, shifts, npts, gwidth, dt)

        f_hat_predicted = _convolve(p, g)

        resid, sumsq_ip1 = get_residual(f_hat, f_hat_predicted)
        sumsq_ip1 = sumsq_ip1 / power
        d_error = 100 * (sumsq_i - sumsq_ip1)

        log.info('%10d  %16.9e  %10.3f   %7.2f%%   %9.4f%%', nshifts, dt * amps[-1], shifts[-1] * dt,
                 100 * sumsq_ip1, d_error)
    # end while

    # *************************************************************************

    log.info('Last Error Change = %9.4f%%', d_error)

    # if the last change made no difference, drop it
    fit = 100 - 100 * sumsq_ip1
    if d_error <= tol:
        nshifts = nshifts - 1
        fit = 100 - 100 * sumsq_i
        log.info('Hit the min improvement tolerance - halting.')
    # end if

    if nshifts >= max_pulses:
        log.warning('Hit the max number of pulses - halting.')
    # end if

    log.info('Number of pulses in final result: %d', nshifts)
    log.info('The final deconvolution reproduces %6.1f%% of the signal.', fit)

    # *************************************************************************

    # compute the final prediction
    p, pulses = build_decon(amps, shifts, npts, gwidth, dt)

    response_trace = denominator.copy()  # clone input trace
    f_hat_predicted = _convolve(p, g)
    response_trace.data = f_hat_predicted

    # Shift RF estimate to synchronize with time line of inputs. If this causes RF to wrap around in the time domain,
    # try increasing the time extent of the input data.
    p_shifted = phase_shift(p, -time_shift, dt)
    rf_trace = numerator.copy()  # clone input trace
    rf_trace.data = p_shifted

    return rf_trace, pulses, f_hat, response_trace, fit

# end func


#------------------------------------------------------------------------------

def main():
    """
    Main
    """
    use_Ammon_data = False
    if use_Ammon_data:
        # Load Ammon test data (has very low sampling rate of 5 Hz)
        src_sac_z = r"C:\software\hiperseis\seismic\receiver_fn\DATA\Ammon_test\lac_sp.z"
        src_sac_r = r"C:\software\hiperseis\seismic\receiver_fn\DATA\Ammon_test\lac_sp.r"
        src_expected_rf = r"C:\software\hiperseis\seismic\receiver_fn\DATA\Ammon_test\lac.i.eqr"
        src_observed = r"C:\software\hiperseis\seismic\receiver_fn\DATA\Ammon_test\observed"
        src_predicted = r"C:\software\hiperseis\seismic\receiver_fn\DATA\Ammon_test\predicted"
        source = obspy.read(src_sac_z, format='sac')
        response = obspy.read(src_sac_r, format='sac')
        expected_rf = obspy.read(src_expected_rf, format='sac')
        observed_radial = obspy.read(src_observed, format='sac')
        predicted_radial = obspy.read(src_predicted, format='sac')
        # Run algorithm on Ammon input
        time_shift = 5.0
        rf_trace, _, _, predicted_response, fit = iter_deconv_pulsetrain(
            source[0], response[0], 200, gwidth=2.5, time_shift=time_shift)

        # Plot the results of expected vs computed
        times = np.arange(len(observed_radial[0]))/observed_radial[0].stats.sampling_rate - time_shift
        plt.figure(figsize=(12, 8))
        plt.plot(times, observed_radial[0].data, alpha=0.8)
        # plt.plot(times, expected_response, ":", alpha=0.8)
        plt.plot(times, predicted_radial[0].data, alpha=0.8)
        plt.plot(times, predicted_response.data, alpha=0.8)
        plt.xlabel("Time (s)")
        plt.ylabel("Radial amplitude")
        plt.text(0.02, 0.02, "Prediction match to observation: {:.2f}%".format(fit), transform=plt.gca().transAxes)
        plt.grid("#80808080", linestyle=':')
        plt.legend(['Observed', 'Predicted observation (Ammon)', 'Predicted observation (Python port)'])
        plt.title("Iterative deconv test, Radial component")
        plt.savefig('ammon_iterdeconvd_observation.png', dpi=300)
        plt.close()

        plt.figure(figsize=(12, 8))
        # Normalize so that sum of squares is 1.
        sum_sq = np.sum(np.square(expected_rf[0].data))
        plt.plot(times, expected_rf[0].data/np.sqrt(sum_sq), '-o', markersize=2, fillstyle='none', alpha=0.8)
        sum_sq = np.sum(np.square(rf_trace.data))
        plt.plot(times, rf_trace.data/np.sqrt(sum_sq), '-^', markersize=2, fillstyle='none', alpha=0.8)
        plt.xlabel("Time (s)")
        plt.ylabel("RF amplitude (arb. units)")
        plt.grid("#80808080", linestyle=':')
        plt.legend(['Ammon reference RF', 'Computed RF (Python port)'])
        plt.title("Iterative deconv test, Receiver Function")
        plt.savefig('ammon_iterdeconvd_RF.png', dpi=300)
        plt.close()

    else:
        # Load test data from Bilby. This is broadband data sampled at 50 Hz, so it is downsampled first before passing
        # to iterative deconvolution.
        src_trace_file = r"C:\software\hiperseis\seismic\receiver_fn\DATA\7W.BL05_event_waveforms_for_rf_filtered.h5"
        src_data = obspy.read(src_trace_file, format='h5')
        # Run deconv on traces associated with same events
        max_iterations = 200
        time_window = (-10, 30)
        i = 0
        all_traces = []
        freq_cutoff = (0.25, 1.0)
        outfolder = 'test_iterdeconv_freq{:.2f}-{:.2f}'.format(*freq_cutoff)
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        # end if
        for stream3c in rf.util.IterMultipleComponents(src_data.copy(), key='onset', number_components=(2, 3)):
            stream3c.filter('bandpass', freqmin=freq_cutoff[0], freqmax=freq_cutoff[1], corners=2,
                            zerophase=True).interpolate(5.0)
            rf_stream = rf.RFStream(stream3c)
            rf_stream.rotate('NE->RT')
            rf_stream.trim2(*time_window, reftime='onset')
            rf_stream.detrend('linear')
            rf_stream.taper(0.2, max_length=0.5)
            source = rf_stream .select(component='Z')[0]
            response = rf_stream .select(component='R')[0]

            rf_trace, pulses, expected_response, predicted_response, fit = iter_deconv_pulsetrain(
                source, response, max_iterations, gwidth=2.5)
            all_traces.append(rf_trace)

            # Generate plots
            event_id = source.stats.event_id
            times = source.times() - (source.stats.onset - source.stats.starttime)
            plt.figure(figsize=(12, 8))
            plt.subplot(211)
            plt.plot(times, expected_response, alpha=0.8)
            plt.plot(times, predicted_response, alpha=0.8)
            plt.xlabel("Time (s)")
            plt.ylabel("Radial amplitude")
            plt.text(0.02, 0.07, "Input filter band: ({:.2f}, {:.2f}) Hz".format(*freq_cutoff),
                     fontsize=8, transform=plt.gca().transAxes, color="#404040")
            plt.text(0.02, 0.02, "Prediction match to observation: {:.2f}%".format(fit), fontsize=8,
                     transform=plt.gca().transAxes, color="#404040")
            plt.grid("#80808080", linestyle=':')
            plt.legend(['Expected R-component', 'Predicted by RF'])
            plt.title("Event {} observed vs predicted Radial component".format(event_id))

            plt.subplot(212)
            sum_sq = np.sum(np.square(rf_trace.data))
            plt.plot(times, rf_trace.data/np.sqrt(sum_sq))
            plt.xlabel("Time (s)")
            plt.ylabel("RF amplitude (arb. units)")
            plt.grid("#80808080", linestyle=':')
            plt.legend(['Computed RF'])
            plt.title("Estimated Receiver Function", y=0.9)

            plt.savefig(os.path.join(outfolder, '{:03d}.png'.format(i)), dpi=300)
            plt.close()

            i += 1
        # end for

        all_rf_stream = rf.RFStream(all_traces).sort(keys=['back_azimuth'])
        stack_file = os.path.join(outfolder, 'rf_stack.png')
        plot_rf_stack(all_rf_stream, time_window=time_window, trace_height=0.12, save_file=stack_file, dpi=300)
        plt.close()
    # end if

# end main


if __name__ == "__main__":
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    main()
