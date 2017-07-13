import datetime
import glob, os
from os.path import join, exists
import obspy
from obspy.core import Stream, UTCDateTime
from obspy import read
from obspy.signal.cross_correlation import xcorr
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
import math
from fft import *

# Greatest common divisor of more than 2 numbers.

def gcd(*numbers):
    """Return the greatest common divisor of the given integers"""
    from fractions import gcd
    return reduce(gcd, numbers)

# Assuming numbers are positive integers.

def lcm(*numbers):
    """Return lowest common multiple."""    
    def lcm(a, b):
        return (a * b) // gcd(a, b)
    return reduce(lcm, numbers, 1)


def butterworthBPF(spectrum, samprate, f_lo, f_hi):
    f = np.fft.fftfreq(spectrum.shape[0], 1.0 / samprate)
    p = 4  # poles
    gain_lpf = 1.0 / (1 + (f / f_hi) ** (2 * p))
    gain_hpf = 1.0 / (1 + (f / f_lo) ** (2 * p))
    gain_bpf = gain_lpf - gain_hpf
    spectrum *= gain_bpf
    return spectrum


def zeropad(tr, padlen):
    assert (tr.shape[0] < padlen)
    padded = np.zeros(padlen)
    padded[0:tr.shape[0]] = tr
    return padded


def zeropad_ba(tr, padlen):
    assert (tr.shape[0] < padlen)
    padded = np.zeros(padlen, dtype=np.complex_)
    s = int((padlen - tr.shape[0]) / 2)
    padded[s:(s + tr.shape[0])] = scipy.fftpack.fftshift(tr)
    return scipy.fftpack.ifftshift(padded)

def zeropad_ba_old(tr, padlen):
    assert (tr.shape[0] < padlen)
    padded = np.zeros(padlen, dtype=np.complex_)
    s = int((padlen - tr.shape[0]) / 2)
    padded[s:(s + tr.shape[0])] = tr
    return padded

def taperDetrend(tr, taperlen):
    """ First detrend, then tapers signal"""
    tr -= np.median(tr)
    tr[0:taperlen] *= 0.5 * (1 + np.cos(np.linspace(-math.pi, 0, taperlen)))
    tr[-taperlen:] *= 0.5 * (1 + np.cos(np.linspace(0, math.pi, taperlen)))
    return tr


def xcorr2(tr1, tr2, window_seconds=3600, interval_seconds=86400, flo=0.9, fhi=1.1):
    # print "About to xcor..."
    sr1 = tr1.stats.sampling_rate
    sr2 = tr2.stats.sampling_rate
    tr1_d_all = tr1.data  # refstn
    tr2_d_all = tr2.data
    lentr1_all = tr1_d_all.shape[0]
    lentr2_all = tr2_d_all.shape[0]
    window_samples_1 = window_seconds * sr1
    window_samples_2 = window_seconds * sr2
    interval_samples_1 = interval_seconds * sr1
    interval_samples_2 = interval_seconds * sr2
    itr1s = 0
    itr2s = 0
    resll = []
    while itr1s < lentr1_all and itr2s < lentr2_all:
        itr1e = min(lentr1_all, itr1s + interval_samples_1)
        itr2e = min(lentr2_all, itr2s + interval_samples_2)
        wtr1s = 0
        wtr2s = 0
        resl = []
        sr = max(sr1, sr2)
	xcorlen = 2 * window_seconds * sr - 1
        fftlen = 2 ** (int(np.log2(xcorlen)) + 1)
        while wtr1s < itr1e and wtr2s < itr2e:
            wtr1e = int(min(itr1e, wtr1s + window_samples_1))
            wtr2e = int(min(itr2e, wtr2s + window_samples_2))
            if wtr1e - wtr1s < window_samples_1 or wtr2e - wtr2s < window_samples_2:
                wtr1s = wtr1e
                wtr2s = wtr2e
                continue
            tr1_d = taperDetrend(tr1_d_all[wtr1s:wtr1e].astype(np.float_, copy=True), int(sr1 / flo))
            tr2_d = taperDetrend(tr2_d_all[wtr2s:wtr2e].astype(np.float_, copy=True), int(sr2 / flo))
            if (sr1 < sr2):
                fftlen2 = fftlen
                fftlen1 = int((fftlen2 * 1.0 * sr1) / sr)
                outdims2 = np.array([fftlen2])
                outdims1 = np.array([fftlen1])
                rf = zeropad_ba(fftn(zeropad(tr1_d, fftlen1), shape=outdims1), fftlen2) * fftn(
                    zeropad(ndflip(tr2_d), fftlen2), shape=outdims2)
            elif (sr1 > sr2):
                fftlen1 = fftlen
                fftlen2 = int((fftlen1 * 1.0 * sr2) / sr)
                outdims2 = np.array([fftlen2])
                outdims1 = np.array([fftlen1])
                rf = fftn(zeropad(tr1_d, fftlen1), shape=outdims1) * zeropad_ba(
                    fftn(zeropad(ndflip(tr2_d), fftlen2), shape=outdims2), fftlen1)
            else:
                fftlen = 2 ** (int(np.log2(2 * window_samples_1 - 1)) + 1)
                outdims = np.array([fftlen])
                rf = fftn(zeropad(tr1_d, fftlen), shape=outdims) * fftn(zeropad(ndflip(tr2_d), fftlen), shape=outdims)
            resl.append(rf)
            wtr1s = wtr1e
            wtr2s = wtr2e
        step = np.sign(np.fft.fftfreq(fftlen, 1.0 / sr))
        mean = reduce((lambda tx, ty: tx + ty), resl) / len(resl)
        mean = butterworthBPF(mean, sr, flo, fhi)
        mean = mean + step * mean  # compute analytic
        mean = ifftn(mean)
        # mask out the zero line for the max
        zerolower = xcorlen / 2 - sr / flo
        zeroupper = xcorlen / 2 + sr / flo
        maskmean = mean.view(np.ma.MaskedArray)
        maskmean[zerolower:zeroupper] = np.ma.masked
        #mean /= np.ma.max(maskmean)
        mean /= 3.0*np.ma.std(maskmean) # three sigma
        mean = np.sqrt(np.imag(mean) ** 2 + np.real(mean) ** 2)
        resll.append(mean[:xcorlen])
        itr1s = wtr1e
        itr2s = wtr2e
    return np.array(resll)


def IntervalStackXCorr(refstn, st,
                       window_seconds=3600,
                       interval_seconds=86400, flo=0.9, fhi=1.1):
    comp_list = []
    a_list = []
    b_list = []
    ref_list = []
    stn_list = []
    xcorr_list = []
    xcorr_x_list = []
    for tr_1 in refstn:
        station_1 = tr_1.stats.station

        for tr_2 in st:
            station_2 = tr_2.stats.station

            comp_list.append(station_2 + '_' + station_1)

            xcorr_func = xcorr2(tr_1, tr_2, window_seconds, interval_seconds, flo, fhi)

            tr_1_len = tr_1.stats.endtime.timestamp - tr_1.stats.starttime.timestamp
            tr_2_len = tr_2.stats.endtime.timestamp - tr_2.stats.starttime.timestamp
            xcorr_len = tr_2_len - tr_1_len
            xcorr_start = tr_2.stats.starttime.timestamp

            X = np.linspace(-window_seconds, window_seconds, xcorr_func.shape[1])
            ref_list.append(tr_1)
            stn_list.append(tr_2)
            xcorr_list.append(xcorr_func)
            xcorr_x_list.append(X)
    return [xcorr_list, xcorr_x_list, comp_list]


def saveXCorr(xcorr_list, xcorr_x_list, xcor_output_dir, figname):
    os.chdir(xcor_output_dir)
    for i, l in enumerate(xcorr_list):
        np.savetxt(figname + '.xcor', xcorr_list[i])


def saveXCorrPlot(xcorr_list, xcorr_x_list, plot_output_dir, figname, comp_list):
    nplots = len(comp_list)
    if nplots == 0:
        return
    fig, ax = plt.subplots(1, nplots, squeeze=False)
    fig.set_figheight(23.39)  # 11.69)
    fig.set_figwidth(16.53)  # 8.27)

    for i, comp_xcorr in enumerate(comp_list):
        ax[0, i].pcolormesh(xcorr_x_list[i], np.arange(xcorr_list[i].shape[0]+1), xcorr_list[i], vmin=0, vmax=1,cmap='jet')  # , label=comp_xcorr, lw=1.2)
        ax[0, i].set_ylabel('Day')
        ax[0, i].set_xlabel('Lag (s)')
        ax[0, i].set_title(comp_xcorr)
        ax[0, i].invert_yaxis()

    os.chdir(plot_output_dir)
    plt.savefig(figname + '.png', dpi=300)
    plt.clf()


if __name__ == "__main__":
    sdr = '/g/data/ha3/Passive/Ref/'
    fn1 = 'AUQLP__BHZ__-14-05-23-00-00-00.msd'
    fn1b = 'AUQLP__BHZ__-14-05-24-00-00-00.msd'
    fn2 = 'AUINKA_BHZ__-14-05-23-00-00-00.msd'
    fn2b = 'AUINKA_BHZ__-14-05-24-00-00-00.msd'
    st1 = read(sdr+fn1)
    st1 += read(sdr+fn1b)
    st1.merge()
    st2 = read(sdr+fn2)
    st2 += read(sdr+fn2b)
    st2.merge()
    st2.resample(80.0,no_filter=True)
    print(st1)
    print(st2)
    ylist, x, comp_list = IntervalStackXCorr(st1,st2)
    print ylist
    print x
    print comp_list
    st1.resample(80.0,no_filter=True)
    ylist2, x2, comp_list2 = IntervalStackXCorr(st1,st2)
    print ylist2
    print x2
    print comp_list2
    # confirm the error is zero
    #stdyratio = np.std(ylist2[0]/ylist[0])
    maxdiff = np.max(np.abs(ylist2[0]-ylist[0]))
    saveXCorrPlot(ylist,x,'/g/data/ha3/','outxcor_dsr4',comp_list)
    saveXCorr(ylist,x2,'/g/data/ha3/','outxcor_8040')
    saveXCorr(ylist2,x2,'/g/data/ha3/','outxcor_8080')
    #print "Std y ratio = " + str(stdyratio)
    print "Max diff = " + str(maxdiff)
    #plt.subplot(211)
    #plt.plot(ylist[0][179999:180001])
    #plt.subplot(212)
    #plt.plot(ylist2[0][179999:180001])
    #plt.show()
    assert maxdiff < 0.05, "The xcor after resampling is not the same"
    
    
