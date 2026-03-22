from collections import defaultdict
import logging
import numpy as np
from scipy.signal import correlate, find_peaks, peak_prominences
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy import stats
import rf
from seismic.receiver_fn import rf_util

# pylint: disable=invalid-name, logging-format-interpolation

logging.basicConfig()

def is_delayed(cha_data, dt_max=0.2):
    """
    Checks if the data shows signs of the presence of reverberations

    :param cha_data: List or iterable of RF traces to use for H-k stacking.
    :type cha_data: Iterable(rf.RFTrace)
    :param dt_max: if the median temporal offset between an RF peak and the onset time > dt_max,
                   this function returns true
    :type dt_max: float

    return: Bool
    """

    dt_array = []
    for i, tr in enumerate(cha_data):
        lead_time = tr.stats.onset - tr.stats.starttime

        relative_time = tr.times() - lead_time
        mask = np.array((relative_time < 0) | (relative_time > 5.0))
        loc = np.argmax(np.ma.masked_array(tr.data, mask=mask))
        dt_array.append(relative_time[loc])
    # end for
    dt_array = np.array(dt_array)

    if(np.median(dt_array) > dt_max): return True

    return False
# end func

class ResonanceFilter:
    def __init__(self, cha_stream, autocorr_window=20):
        """
        @param cha_stream: RFStream containing the R/L channel
        @param autocorr_window: window (seconds) over which autocorrelations are computed
        Creates a resonance filter and associated attributes
        """
        self.autocorr_window = autocorr_window

        stackable_stream = rf_util.get_stackable_stream(cha_stream)

        num_stackable = len(stackable_stream)
        if num_stackable < len(cha_stream):
            print('Removed {} traces from in ResonanceFilter to make stream stackable!'.format(num_stackable))
        # end if

        # compute mean RF
        mtr = stackable_stream[0].copy()
        mtrd = np.zeros(mtr.data.shape)
        for tr in stackable_stream: mtrd += tr.data
        mtrd /= len(stackable_stream)
        mtr.data = mtrd

        tr = mtr.copy()
        tr = tr.slice(tr.stats.onset, tr.stats.onset + autocorr_window)
        times = tr.times()

        # compute autocorrelation
        af = correlate(tr.data, tr.data, mode='full')
        af = af[len(af) // 2:]

        afn = af / np.max(af) # normalize autocorrelation
        r0 = -(np.min(afn))
        Dt = np.argmin(afn) * 1. / tr.stats.sampling_rate

        # assemble resonance filter:
        resonance_filter = np.zeros(len(mtr.data))
        resonance_filter[0] = 1
        resonance_filter[int(Dt * mtr.stats.sampling_rate)] = r0

        # generate modelled af
        def func(t, c, a, p):
            return c * np.exp(-a * t) * np.cos(np.pi * t / p)

        # end func
        opt, _ = curve_fit(func, times, af)
        af_fitted = func(times, *opt)

        if(0):
            fig, ax = plt.subplots()
            fig.set_size_inches(20, 10)
            ax.plot(times, af, c='k', lw='2')
            ax.plot(times, func(times, *opt), c='r')
            print('r0: {}, Dt: {}'.format(r0, Dt))
        # end if

        # save attributes
        self.r0 = r0
        self.Dt = Dt
        self.filter = resonance_filter
        self.mean_trace = mtr
        self.fitted_mean_af = af_fitted
    # end func
# end class

def has_reverberations(cha_stream, autocorr_window=20):
    """
    Checks if reverberations are present, based on criteria outlined in
    Cunningham and Lekic 2019.
    @param cha_stream: RF channel stream
    @param autocorr_window: window (seconds) over which autocorrelations are
           performed while designing resonance filter
    @return: boolean
    """
    if not is_delayed(cha_stream): return False

    rfilter = ResonanceFilter(cha_stream, autocorr_window=autocorr_window)
    mtr = rfilter.mean_trace

    # find the two most prominent peaks in descending order
    peaks, _ = find_peaks(mtr.data, height=0)
    prominences = peak_prominences(mtr.data, peaks)[0]
    peaks = peaks[np.argsort(prominences)[::-1]]

    if(0):
        fig, ax = plt.subplots()
        ax.plot(mtr.times(), mtr)
        ax.plot(mtr.times()[peaks[:2]], mtr.data[peaks[:2]], "x")
        ax.set_xlim(40, 70)
    # end if

    if (len(peaks) >= 2):
        # assess criteria 2 and 3 in Cunningham and Lekic 2019
        C2 = mtr.data[peaks[1]] > mtr.data[peaks[0]] * 0.2  # Ps is at least 20% of the P arrival amplitude
        C3 = mtr.data[peaks[1]] > mtr.data[peaks[0]] * 0.9

        #print('C3: {}'.format(C3))
        if (C3): return True
    else:
        return False
    # end if

    in_stream = []
    out_stream = []
    for tr in cha_stream:
        tr_copy = tr.copy().slice(tr.stats.onset, tr.stats.onset + autocorr_window)
        in_stream.append(tr_copy.copy())

        tr_copy.data = np.convolve(tr_copy.data, rfilter.filter[0:len(tr_copy.data)], mode='full')
        tr_copy.data = tr_copy.data[:len(tr_copy.data) // 2 + 1]

        out_stream.append(tr_copy)
    # end for

    diffs1 = []
    diffs2 = []
    af_fitted = rfilter.fitted_mean_af / np.max(rfilter.fitted_mean_af)
    for itr, otr in zip(in_stream, out_stream):
        if (len(itr) == len(otr)):
            # differences between input and dereverberated traces
            diffs1.append(itr.data / np.max(itr.data) - otr.data / np.max(otr.data))

            # compute trace autocorrelation
            af = correlate(itr.data, itr.data, mode='full')
            af = af[len(af) // 2:]

            # differences between trace autocorrelation and fitted
            # autocorrelation of mean trace
            diffs2.append(af / np.max(af) - af_fitted)
        # end if
    # end for
    diffs1 = np.array(diffs1)
    diffs2 = np.array(diffs2)

    # computre variances of mean differences
    # across all input traces
    v1 = np.var(np.mean(diffs1, axis=0))
    v2 = np.var(np.mean(diffs2, axis=0))

    # assess criteria 1 in Cunningham and Lekic 2019
    C1 = v1 > v2

    #print('C1: {}, C2:{}'.format(C1, C2))
    if (C1 and C2):
        return True
    else:
        return False
# end func

def apply_reverberation_filter(cha_stream, autocorr_window=20):
    """
    Applies dereverberation filter to input traces inplace
    @param cha_stream: input RFsStream
    @param autocorr_window: window (seconds) over which autocorrelations are computed
    @return:
    """
    rfilter = ResonanceFilter(cha_stream, autocorr_window=autocorr_window)

    # apply resonance filter to all traces regardless of length
    result_stream = []
    for tr in cha_stream:
        tr_copy = tr.copy()

        temp = np.convolve(tr_copy.data, rfilter.filter, mode='full')
        tr_copy.data = temp[:len(tr_copy.data)]

        if(tr.data.shape != tr_copy.data.shape):
            print(tr.data.shape, tr_copy.data.shape, tr.id)
        assert tr.data.shape == tr_copy.data.shape, 'Input/output length mismatch detected in ' \
                                                    'reverberation removal routine'

        # find dt
        lead_time = tr.stats.onset - tr.stats.starttime
        relative_time = tr.times() - lead_time
        mask = np.array((relative_time < 0) | (relative_time > 5.0))
        loc = np.argmax(np.ma.masked_array(tr.data, mask=mask))
        dt = relative_time[loc]

        tr_copy.stats.update({'t1_offset': dt,
                              't2_offset': rfilter.Dt - dt,
                              't3_offset': rfilter.Dt})

        result_stream.append(tr_copy)
    # end for

    return rf.RFStream(result_stream)
# end func
