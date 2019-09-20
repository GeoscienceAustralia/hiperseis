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

    # n2 = 1
    # while n2 < n:
    #     n2 = n2 * 2

    # Use the Numerical Recipes routine to compute the cross correlation
    #     call correl(f, g, n2, c) # output ends up in c
    c = scipy.signal.correlate(f, g)
    # Zero or positive lags in result (using sign convention of correl)
    # are in elements [n/2:n] of c.
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

    # n2 = 1
    # while n2 < n:
    #     n2 = n2 * 2
    #
    # halfpts = n2 // 2

    #     call realft(x,halfpts,forward)
    fft_x = np.fft.rfft(x)  # complex array
    n2 = len(fft_x)

    df = 1.0 / (float(n) * dt)
    d_omega = two_pi * df
    gwidth = 4.0 * gwidth_factor * gwidth_factor

    omega = np.arange(n2) * d_omega
    gauss = np.exp(-omega * omega / gwidth)
    fft_x = fft_x * gauss

    #     call realft(x,halfpts,inverse)
    x_filt = np.fft.irfft(fft_x, n)  # real_array

    # scalefactor = dt * (2 * df)
    # x_filt = x_filt * scalefactor

    return x_filt


def phs_shift(x, theshift, dt):
    """
    phase shifts a signal
    """
    n = len(x)
    two_pi = 2 * np.pi

    # n2 = 1
    # while n2 < n:
    #     n2 = n2 * 2

    # halfpts = n2 // 2

#     call realft(x, halfpts, forward)
    fft_x = np.fft.rfft(x)  # complex array
    n2 = len(fft_x)

    df = 1 / (float(n) * dt)
    d_omega = two_pi * df

    omega = np.arange(n2) * d_omega
    # The phase shift is omega (angular velocity) * delta_t (the time shift).
    spectral_shift = np.exp(omega*theshift*1j)
    fft_x = fft_x * spectral_shift

#     call realft(x,halfpts,inverse)
    x_shifted = np.fft.irfft(fft_x, n)  # real_array

    # scalefactor = dt * (2 * df)
    # x_shifted = x_shifted * scalefactor

    return x_shifted


def build_decon(amps, shifts, n, gwidth, dt):
    """
    compute the predicted time series from a set of
    amplitudes and shifts
    """
    p = np.zeros(n)
    for i, amp in zip(shifts, amps):
        p[i] += amp
    # p = np.roll(p, n//2)

    p_filt = gauss_filter(p, gwidth, dt)

    return p_filt


def _convolve(x, y, dt):
    """
    """
    n = len(x)
    assert len(y) == n
    z = scipy.signal.convolve(x, y)
    z_sync = z[:n]
    return z_sync
    # z = scipy.signal.convolve(x, y)*2*dt
    # z_pos = z[-n:]
    # return z_pos


def iterdeconfd(denominator, numerator, maxbumps, tol=1.0e-3, gwidth=2.0, time_shift=None, lpositive=False,
                verbose=False):
    """
    Iterative deconvolution of source and response signal to generate seismic receiver function.

    :param denominator: The source signal
    :param response: The response signal (e.g. R or T component)
    :return:
    """
    MAXPTS = 100000
    assert len(denominator) <= MAXPTS  # Sanity check input size
    assert len(numerator) <= MAXPTS
    MAXG = 200  # Maximum number of pulses to synthesize in p

    amps = []
    shifts = []

    print()
    print('Program iterdeconfd - Version 1.04, 1997-98')
    print('Chuck Ammon, Saint Louis University')
    print('Adapted to Python by Andrew Medlin, Geoscience Australia (2019)')
    print()

    if maxbumps > MAXG:
        print('Maximum Number of bumps is %d' % MAXG)
        maxbumps = MAXG
    # end if

    # theshift = float(input('What is the phase shift (secs) for the output?'))
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
    # nptsd = len(g)
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
    # xc.write('ccor0.sac', format='sac')

    # find the peak in the correlation
    # maxlag = npts // 2
    # print('The maximum spike delay is %g sec' % (float(maxlag) * dt))

    # g_prime = g[0:maxlag]
    if lpositive:
        shifts.append(np.argmax(xc))
    else:
        xc_abs = np.abs(xc)
        shifts.append(np.argmax(xc_abs))
    # end if
    amps.append(xc[shifts[-1]])

    nshifts = 0

    # compute the predicted deconvolution result
    p = build_decon(amps, shifts, npts, gwidth, dt)
    if verbose:
        p_shifted = phs_shift(p, -time_shift, dt)
        p_shifted.write('d000.sac', format='sac')
    # end if

    # convolve the prediction with the denominator signal
    f_hat_predicted = _convolve(p, g, dt)

    if verbose:
        f_hat_predicted.write('p000.sac', format='sac')
    # end if

    if verbose:
        resfile = 'r%03d.sac' % 0
        f_hat.write(resfile, format='sac')
    # end if

    # compute the residual (initial error is 1.0)
    resid, sumsq_ip1 = get_residual(f_hat, f_hat_predicted)

    sumsq_i = 1.0
    sumsq_ip1 = sumsq_ip1 / power
    d_error = 100 * (sumsq_i - sumsq_ip1)

    resfile = 'r%03d.sac' % 0
    if verbose:
        resid.write(resfile, format='sac')
    # end if

    print('%11s %s' % ('File', '  Spike amplitude  Spike delay   Misfit   Improvement'))
    print('%10s  %16.9e  %10.3f   %7.2f%%   %9.4f%%' % (
        resfile, dt * amps[-1], (shifts[-1] - 1) * dt, 100 * sumsq_ip1, d_error))

    # *************************************************************************
    while (np.abs(d_error) > tol) and (nshifts < maxbumps):

        nshifts = nshifts + 1
        sumsq_i = sumsq_ip1

        xc = fcorrelate(resid, g_hat)

        # g_prime = g[0:maxlag]
        if lpositive:
            shifts.append(np.argmax(xc))
        else:
            xc_abs = np.abs(xc)
            shifts.append(np.argmax(xc_abs))
        # end if
        amps.append(xc[shifts[-1]])

        p = build_decon(amps, shifts, npts, gwidth, dt)
        if verbose:
            filename = 'd%03d.sac' % nshifts
            p_shifted = phs_shift(p, -time_shift, dt)
            # wsac1(filename, p, npts, -theshift, dt, nerr)
            p_shifted.write(filename, format='sac')
        # end if

        f_hat_predicted = _convolve(p, g, dt)
        if verbose:
            filename = 'p%03d.sac' % nshifts
            f_hat_predicted.write(filename, format='sac')
        # end if

        resid, sumsq_ip1 = get_residual(f_hat, f_hat_predicted)

        sumsq_ip1 = sumsq_ip1 / power
        resfile = 'r%03d.sac' % nshifts
        if verbose:
            resid.write(resfile, format='sac')
        # end if
        d_error = 100 * (sumsq_i - sumsq_ip1)

        print('%10s  %16.9e  %10.3f   %7.2f%%   %9.4f%%' % (resfile,
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
    p = build_decon(amps, shifts, npts, gwidth, dt)

    trace = denominator.copy()  # clone input trace
    f_hat_predicted = _convolve(p, g, dt)
    trace.data = f_hat_predicted
    trace.write('predicted.sac', format='sac')

    # write out the answer
    # p = build_decon(amps, shifts, npts, gwidth, dt)
    p_shifted = phs_shift(p, -time_shift, dt)

    # This rigmarole is to do with writing detailed SAC header information to file
    #       call newhdr
    #       call rsac1(numerator, g, ndummy, b, dt, MAXPTS, nerr)
    #       call setnhv('NPTS',npts,nerr)
    #       call setfhv('B',-theshift,nerr)
    #       theend = -thshift + (npts-1)*dt
    #       call setfhv('E',theend,nerr)
    #       call setnhv('NZSEC',-12345,nerr)
    #       call setfhv('USER0',gwidth,nerr)
    # c     call setkhv('KUSER0','Rftn',nerr)
    # c     call setkhv('KUSER1','IT_DECON',nerr)
    #       call wsac0('decon.out',xdummy,p,nerr)
    trace = numerator.copy()  # clone input trace
    trace.data = p_shifted
    trace.write('decon.sac', format='sac')

    # write out the gaussian filter
    if verbose:
        #         call newhdr
        #         call zero(p,MAXPTS)
        #         p(1) = 1 / dt
        #         call phs_shift(p,theshift,npts,dt)
        #         call gfilter(p,gwidth,npts,dt)
        #         call wsac1('thefilter',p,npts,beg,dt,nerr)
        p = np.zeros(npts)
        p[0] = 1.0 / dt
        p_shifted = phs_shift(p, -time_shift, dt)
        p = gauss_filter(p_shifted, gwidth, dt)
        p.write('thefilter.sac', format='sac')
    # end if

# end func

def iter3c(stream):
    return rf.util.IterMultipleComponents(stream, key='onset', number_components=(2, 3))

#------------------------------------------------------------------------------

def main():
    """
    Main
    """
    src_sac_z = r"C:\software\hiperseis\seismic\receiver_fn\DATA\Ammon_test\lac_sp.z"
    src_sac_r = r"C:\software\hiperseis\seismic\receiver_fn\DATA\Ammon_test\lac_sp.r"
    src_expected_rf = r"C:\software\hiperseis\seismic\receiver_fn\DATA\Ammon_test\lac.i.eqr"
    source = obspy.read(src_sac_z, format='sac')
    response = obspy.read(src_sac_r, format='sac')
    expected_rf = obspy.read(src_expected_rf, format='sac')
    result = iterdeconfd(source[0], response[0], 200, gwidth=2.5, time_shift=8.0)

    # src_trace_file = r"C:\software\hiperseis\seismic\receiver_fn\DATA\7W.BL05_event_waveforms_for_rf_filtered.h5"
    # src_data = obspy.read(src_trace_file, format='h5')
    # max_iterations = 200
    # time_window = (-10, 30)
    # for stream3c in iter3c(src_data.copy()):
    #     rf_stream = rf.RFStream(stream3c)
    #     rf_stream .rotate('NE->RT')
    #     rf_stream .trim2(*time_window, reftime='onset')
    #     rf_stream.detrend('linear')
    #     rf_stream.taper(0.2, max_length=0.5)
    #     source = rf_stream .select(component='Z')[0]
    #     response = rf_stream .select(component='R')[0]
    #     result = iterdeconfd(source, response, max_iterations, gwidth=10.0)

# end main

if __name__ == "__main__":
    main()