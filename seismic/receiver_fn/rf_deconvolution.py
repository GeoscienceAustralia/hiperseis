#!/usr/bin/env python

import numpy as np
import scipy
import scipy.signal

import obspy
import rf

from seismic.receiver_fn.rf_h5_file_event_iterator import IterRfH5FileEvents

def fcorrelate(f, g, n, MAXPTS, dt):  # n and MAXPTS probably redundant
    """
    correlation routine - correlates f and g and replaces the
      g with the cross-correlation the value is normalized
      by the zero-lag autocorrelation of g

    n = len(g), should be same as len(f)

    Returns result in g
    """
    #     real f(MAXPTS), g(MAXPTS), c(8192)
    c = np.zeros(MAXPTS)

    # compute the zero-lag autocorrelation of g
    sum0 = dt * np.dot(g, g)

    n2 = 1
    while n2 < n:
        n2 = n2 * 2

    # Use the Numerical Recipes routine to compute the cross correlation
    #     call correl(f, g, n2, c) # output ends up in c
    c = scipy.signal.correlate(f, g, mode='same')
    # Zero or positive lags in result (using sign convention of correl)
    # are in elements [n/2:n] of c. Roll result to match layout out output
    # of correl.
    c = np.roll(c, -n // 2)

    temp = dt / sum0
    # **Note**: this is intended to change values in g in-place up to n2,
    # i.e. re-using g as output array.
    g[0:n2] = c[0:n2] * temp


def getres(x, y, n):
    """
    get the residual between x and y

    n = len(x), len(y)
    """
    #     real x(n), y(n), r(n), sumsq

    r = x - y
    sumsq = np.dot(r, r)

    return r, sumsq


def gfilter(x, gwidth_factor, n, dt):
    """
    convolve a function with a unit-area Gaussian filter.
    """
    #     real x(n), two_pi, gauss, d_omega, omega
    #     real gwidth, gwidth_factor, sum

    two_pi = 2 * np.pi

    # n2 = 1
    # while n2 < n:
    #     n2 = n2 * 2
    #
    # halfpts = n2 // 2

    #     call realft(x,halfpts,forward)
    fft_x = np.fft.rfft(x)  # complex array
    n2 = len(fft_x)

    df = 1.0 / (float(n2) * dt)
    d_omega = two_pi * df
    gwidth = 4.0 * gwidth_factor * gwidth_factor

    omega = np.arange(n2) * d_omega
    gauss = np.exp(-omega * omega / gwidth)
    fft_x = fft_x * gauss

    #     call realft(x,halfpts,inverse)
    x_filt = np.fft.irfft(fft_x)  # real_array

    scalefactor = dt * (2 * df)
    x_filt = x_filt * scalefactor

    return x_filt


def phs_shift(x, theshift, n, dt):
    """
    phase shifts a signal
    """
#     real x(n), pi, two_pi, theshift, d_omega, omega

    two_pi = 2 * np.pi

    n2 = 1
    while n2 < n:
        n2 = n2 * 2

    halfpts = n2 // 2

#     call realft(x, halfpts, forward)
    fft_x = np.fft.rfft(x, halfpts)  # complex array

    df = 1 / (float(n2) * dt)
    d_omega = two_pi * df

    omega = np.arange(halfpts)*d_omega
    # Copied from Ammon, but weird(?).
    # Shouldn't it just be spectral_shift = np.exp(theshift*1j) ?? (drop the omega)
    spectral_shift = np.exp(omega*theshift*1j)
    ffx_x = fft_x * spectral_shift

#     call realft(x,halfpts,inverse)
    x_shifted = np.fft.irfft(fft_x, halfpts)  # real_array

    scalefactor = dt * (2 * df)
    x_shifted = x_shifted * scalefactor

    return x_shifted


def build_decon(amps, shifts, nshifts, p, n, gwidth, dt):
    """
    compute the predicted time series from a set of
    amplitudes and shifts
    """
    #     real p(n), amps(nshifts)
    #     integer shifts(nshifts)
    #     integer i, n, nshifts

    p[shifts] = p[shifts] + amps

    x = gfilter(p, gwidth, n, dt)

    return x


def iterdeconfd(denominator, numerator, maxbumps, tol=1.0e-3, gwidth=2.0, lpositive=False, verbose=False):
    """
    Iterative deconvolution of source and response signal to generate seismic receiver function.

    :param denominator: The source signal
    :param response: The response signal (e.g. R or T component)
    :return:
    """
    MAXPTS = 8192  # TODO: Update or remove this limit
    assert len(denominator) <= MAXPTS
    assert len(numerator) <= MAXPTS
    MAXG = 200  # TODO: Relax this limit
    f = np.zeros(MAXPTS)
    g = np.zeros(MAXPTS)
    p = np.zeros(MAXPTS)
    #     r = np.zeros(MAXPTS)
    amps = np.zeros(MAXG)
    shifts = np.zeros(MAXG).astype(int)
    resfile = ' ' * 12
    filename = ' ' * 12

    # stdin = 5
    # stdout = 6
    # ounit = 9
    # inunit = 10
    #     forward = 1
    #     inverse = -1
    # lpositive = False
    # verbose = False
    # gwidth = 2.5

    print('\n')
    print('Program iterdeconfd - Version 1.04, 1997-98')
    print('Chuck Ammon, Saint Louis University')
    print('Adapted to Python by Andrew Medlin, Geoscience Australia (2019)')
    print('\n')

    # read in the names of the input files

    # numerator = input('What is the numerator file?')
    # denominator = input('What is the denominator file?')
    # maxbumps = int(input('What is the max number of iterations?'))

    if (maxbumps > MAXG):
        print('Maximum Number of bumps is %d' % MAXG)
        maxbumps = MAXG
    # end if

    # theshift = float(input('What is the phase shift (secs) for the output?'))
    theshift = float(denominator.stats.onset - denominator.stats.starttime)
    # tol = float(input('What is minimum percent error increase to accept?'))
    # gwidth = float(input('What is is the Gaussian filter width factor?'))
    # idum = int(input('Allow negative pulses? (1->y, 0->no)'))
    #
    # if (idum == 1):
    #     lpositive = False
    # else:
    #     lpositive = True
    # # end if

    # idum = int(input('Minimal (0) or verbose output(1)?'))
    # if (idum == 0):
    #     verbose = False
    # else:
    #     verbose = True
    # # end if

    # *************************************************************************
    #     call rsac1(numerator, f, npts, beg, delta, MAXPTS, nerr)
    f = numerator.copy().data
    npts = len(f)

    #     call rsac1(denominator,g,nptsd,b,dt,MAXPTS,nerr)
    g = denominator.copy().data
    # nptsd = len(g)
    dt = 1.0 / denominator.meta.sampling_rate

    # *************************************************************************
    # Find the next power of two greater than the data
    #  dimensions - use the numerator, zero pad
    n = 1
    while n < npts:
        n = n * 2

    if (n > MAXPTS):
        print('Too many points needed.')
        print('n = %d' % n)
        exit(1)
    # end if

    # *************************************************************************
    # zero-pad the data
    npts = n

    # *************************************************************************
    # FINISHED READING THE FILES
    # *************************************************************************

    # Now begin the cross-correlation procedure

    # Put the filter in the signals
    f = gfilter(f, gwidth, npts, dt)
    g = gfilter(g, gwidth, npts, dt)
    # Buffering to disk for later read-back?
    # f.write('numerator', format='sac')
    # f.write('observed', format='sac')
    # g.write('denominator', format='sac')

    # compute the power in the "numerator" for error scaling
    power = np.dot(f, f)

    # correlate the signals
    fcorrelate(f, g, npts, MAXPTS, dt)
    #     g.write('ccor0', format='sac')

    # find the peak in the correlation
    maxlag = npts // 2
    print('The maximum spike delay is %g' % float(maxlag) * dt)

    g_prime = g[0:maxlag]
    if (lpositive):
        shifts[0] = np.argmax(g_prime)
        amps[0] = g_prime[shifts[0]]
    else:
        g_abs = np.abs(g_prime)
        shifts[0] = np.argmax(g_abs)
        amps[0] = g_abs[shifts[0]]
    # end if
    amps[0] = amps[0] / dt

    nshifts = 0

    # read in the signals again  # TODO: rationalize the I/O in this routine
    #     call rsac1(numerator,f,ndummy,beg,delta,MAXPTS,nerr)
    f = numerator.copy().data
    #     call rsac1(denominator,g,ndummy,b,dt,MAXPTS,nerr)
    g = denominator.copy().data
    dt = 1.0 / g[0].meta.sampling_rate  # grab dt from first trace

    # compute the predicted deconvolution result
    p = np.zeros(MAXPTS)
    p = build_decon(amps, shifts, nshifts, p, npts, gwidth, dt)
    if (verbose):
        p_shifted = phs_shift(p, theshift, npts, dt)
        p_shifted.write('d000', format='sac')
    # end if

    # convolve the prediction with the denominator signal
    #     call convolve(p,g,npts,dt)
    p = scipy.signal.convolve(p, g, mode='same')  # need to check dt scaling here...

    if (verbose):
        p.write('p000', format='sac')
    # end if

    # filter the signals
    f = gfilter(f, gwidth, npts, dt)
    g = gfilter(g, gwidth, npts, dt)

    if (verbose):
        resfile = 'r%03d' % 0
        f.write(resfile, format='sac')
    # end if

    # compute the residual (initial error is 1.0)
    r, sumsq_ip1 = getres(f, p, npts)

    sumsq_i = 1.0
    sumsq_ip1 = sumsq_ip1 / power
    d_error = 100 * (sumsq_i - sumsq_ip1)

    resfile = 'r%03d' % 1
    if (verbose):
        r.write(resfile, format='sac')
    # end if

    print('%9s %s' % ('File', 'Spike amplitude   Spike delay   Misfit   Improvement'))
    print('%10s  %16.9e  %10.3f   %7.2f%%   %9.4f%%' % (
    resfile, dt * amps[0], (shifts[0] - 1) * dt, 100 * sumsq_ip1, d_error))

    # *************************************************************************
    while (d_error > tol) and (nshifts < maxbumps):

        nshifts = nshifts + 1
        sumsq_i = sumsq_ip1

        #         g = np.zeros(MAXPTS)
        g = obspy.read(denominator, format='sac')

        g = gfilter(g, gwidth, npts, dt)
        fcorrelate(r, g, npts, MAXPTS, dt)

        g_prime = g[0:maxlag]
        if (lpositive):
            shifts[nshifts] = np.argmax(g_prime)
            amps[nshifts] = g_prime[shifts[nshifts]]
        else:
            g_abs = np.abs(g_prime)
            shifts[nshifts] = np.argmax(g_abs)
            amps[nshifts] = g_abs[shifts[nshifts]]
        # end if
        amps[nshifts] = amps[nshifts] / dt

        p = np.zeros(MAXPTS)
        p = build_decon(amps, shifts, nshifts, p, npts, gwidth, dt)
        if (verbose):
            filename = '(d%03d)' % nshifts
            p_shifted = phs_shift(p, theshift, npts, dt)
            #             wsac1(filename, p, npts, -theshift, dt, nerr)
            p_shifted.write(filename, format='sac')
        # end if

        #         g = np.zeros(MAXPTS)
        g = obspy.read(denominator, format='sac')
        p = scipy.signal.convolve(p, g, mode='same')  # need to check dt scaling here...
        if (verbose):
            filename = '(p%03d)' % nshifts
            p.write(filename, format='sac')
        # end if

        #         f = np.zeros(MAXPTS)
        f = obspy.read(numerator)
        f = gfilter(f, gwidth, npts, dt)
        r, sumsq_ip1 = getres(f, p, npts)

        sumsq_ip1 = sumsq_ip1 / power
        resfile = '(r%03d)' % nshifts
        if (verbose):
            r.write(resfile, format='sac')
        # end if
        d_error = 100 * (sumsq_i - sumsq_ip1)

        print('%10s  %16.9e  %10.3f   %7.2f%%   %9.4f%%' % (resfile,
                                                            dt * amps[nshifts], (shifts[nshifts] - 1) * dt,
                                                            100 * sumsq_ip1, d_error))
    # enddo

    # *************************************************************************

    print('Last Error Change = %9.4f%%' % d_error)

    # if the last change made no difference, drop it
    fit = 100 - 100 * sumsq_ip1
    if (d_error <= tol):
        nshifts = nshifts - 1
        fit = 100 - 100 * sumsq_i
        print('Hit the min improvement tolerance - halting.')
    # end if

    if (nshifts >= maxbumps):
        print('Hit the max number of bumps - halting.')
    # end if

    print('Number of bumps in final result: %d' % nshifts)
    print('The final deconvolution reproduces %6.1f%% of the signal.' % fit)

    # *************************************************************************

    # compute the final prediction
    p = np.zeros(MAXPTS)
    p = build_decon(amps, shifts, nshifts, p, npts, gwidth, dt)

    #     g = np.zeros(MAXPTS)
    g = obspy.read(denominator, format='sac')
    p = scipy.signal.convolve(p, g, mode='same')  # need to check dt scaling here...
    p.write('predicted', format='sac')
    g = np.zeros(MAXPTS)

    # write out the answer
    p = np.zeros(MAXPTS)
    p = build_decon(amps, shifts, nshifts, p, npts, gwidth, dt)
    p_shifted = phs_shift(p, theshift, npts, dt)

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
    p_shifted.write('decon.out', format='sac')

    # write out the gaussian filter
    if (verbose):
        #         call newhdr
        #         call zero(p,MAXPTS)
        #         p(1) = 1 / dt
        #         call phs_shift(p,theshift,npts,dt)
        #         call gfilter(p,gwidth,npts,dt)
        #         call wsac1('thefilter',p,npts,beg,dt,nerr)
        p = np.zeros(MAXPTS)
        p[0] = 1.0 / dt
        p_shifted = phs_shift(p, theshift, npts, dt)
        p = gfilter(p_shifted, gwidth, npts, dt)
        p.write('thefilter', format='sac')
    # end if

# end func

def iter3c(stream):
    return rf.util.IterMultipleComponents(stream, key='onset', number_components=(2, 3))

#------------------------------------------------------------------------------

def main():
    """
    Main
    """
    src_trace_file = "/g/data/ha3/am7399/rf_validation/7W_bilby/7W.BL05_event_waveforms_for_rf_filtered.h5"
    src_data = obspy.read(src_trace_file, format='h5')
    max_iterations = 200
    time_window = (-10, 30)
    for stream3c in iter3c(src_data.copy()):
        rf_stream = rf.RFStream(stream3c)
        rf_stream .rotate('NE->RT')
        rf_stream .trim2(*time_window, reftime='onset')
        rf_stream.detrend('linear')
        source = rf_stream .select(component='Z')[0]
        response = rf_stream .select(component='R')[0]
        result = iterdeconfd(source, response, max_iterations)
# end main

if __name__ == "__main__":
    main()
