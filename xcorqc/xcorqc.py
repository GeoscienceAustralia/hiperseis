import datetime
import glob, os
from os.path import join, exists
import obspy
from obspy.core import Stream, UTCDateTime
from obspy import read
from obspy.signal.cross_correlation import xcorr
import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
from fft import *


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
                fftlen2 = 2 ** (int(np.log2(2 * max(window_samples_1, window_samples_2) - 1)) + 1)
                fftlen1 = int((fftlen2 * 1.0 * sr1) / sr2)
                fftlen = fftlen2
                outdims2 = np.array([fftlen2])
                outdims1 = np.array([fftlen1])
                rf = zeropad_ba(fftn(zeropad(tr1_d, fftlen1), shape=outdims1), fftlen2) * fftn(
                    zeropad(ndflip(tr2_d), fftlen2), shape=outdims2)
                #rf = zeropad_ba(fftn(zeropad(tr1_d, fftlen1), shape=outdims1), fftlen2) * np.conjugate(fftn(
                #    zeropad(tr2_d, fftlen2), shape=outdims2))
            elif (sr1 > sr2):
                fftlen1 = 2 ** (int(np.log2(2 * max(window_samples_1, window_samples_2) - 1)) + 1)
                fftlen2 = int((fftlen1 * 1.0 * sr2) / sr1)
                fftlen = fftlen1
                outdims2 = np.array([fftlen2])
                outdims1 = np.array([fftlen1])
                rf = fftn(zeropad(tr1_d, fftlen1), shape=outdims1) * zeropad_ba(
                    fftn(zeropad(ndflip(tr2_d), fftlen2), shape=outdims2), fftlen1)
                #rf = fftn(zeropad(tr1_d, fftlen1), shape=outdims1) * np.conjugate(zeropad_ba(
                #    fftn(zeropad(tr2_d, fftlen2), shape=outdims2), fftlen1))
            else:
                fftlen = 2 ** (int(np.log2(2 * window_samples_1 - 1)) + 1)
                outdims = np.array([fftlen])
                rf = fftn(zeropad(tr1_d, fftlen), shape=outdims) * fftn(zeropad(ndflip(tr2_d), fftlen), shape=outdims)
                #rf = fftn(zeropad(tr1_d, fftlen), shape=outdims) * np.conjugate(fftn(zeropad(tr2_d, fftlen), shape=outdims))
            resl.append(rf)
            wtr1s = wtr1e
            wtr2s = wtr2e
        xcorlen = 2 * max(window_samples_1, window_samples_2) - 1
        fftlen = 2 ** (int(np.log2(xcorlen)) + 1)
        sr = max(sr1, sr2)
        step = np.sign(np.fft.fftfreq(fftlen, 1.0 / sr))
        mean = reduce((lambda tx, ty: tx + ty), resl) / len(resl)
        mean = butterworthBPF(mean, tr1.stats.sampling_rate, flo, fhi)
        mean = mean + step * mean  # compute analytic
        mean = ifftn(mean)
        #mean = np.sqrt(np.imag(mean) ** 2 + np.real(mean) ** 2)
        # mask out the zero line for the max
        zerolower = xcorlen / 2 - sr * 1.5
        zeroupper = xcorlen / 2 + sr * 1.5
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
    # for i, comp_xcorr in enumerate(comp_list):
    for i, l in enumerate(xcorr_list):
        np.savetxt(figname + '.xcor', xcorr_list[i])
    # np.savetxt('%s_%d.xcor'%(comp_xcorr,i),xcorr_list[i])


def saveXCorrPlot(xcorr_list, xcorr_x_list, plot_output_dir, figname, comp_list):
    nplots = len(comp_list)
    # print "Number of plots: " + str(nplots)
    if nplots == 0:
        return
    # plt.rc('figure', figsize=(11.69,8.27))
    fig, ax = plt.subplots(1, nplots, squeeze=False)
    fig.set_figheight(23.39)  # 11.69)
    fig.set_figwidth(16.53)  # 8.27)
    # axes = fig.axes
    # print ax

    for i, comp_xcorr in enumerate(comp_list):
        # print xcorr_list[i]
        # ax[i,0].plot(xcorr_x_list[i], (xcorr_list[i]), label=comp_xcorr, lw=1.2)
        # print xcorr_list[i].shape
        ax[0, i].pcolormesh(xcorr_x_list[i], np.arange(xcorr_list[i].shape[0]+1), xcorr_list[i], vmin=0, vmax=1,cmap='jet')  # , label=comp_xcorr, lw=1.2)
        ax[0, i].set_ylabel('Day')
        ax[0, i].set_xlabel('Lag (s)')
        ax[0, i].set_title(comp_xcorr)
        # ax[0,i].grid()
        ax[0, i].invert_yaxis()
    # ax[0,i].legend()
    # ax[1,i].plot(getStreamX(ref_list[i]), (ref_list[i]), label="Reference", lw=1.2)
    # ax[2,i].plot(getStreamX(stn_list[i]), (stn_list[i]), label="Station", lw=1.2)

    # plt.legend()
    # plt.show()
    os.chdir(plot_output_dir)
    # print 'Saving figure...'
    plt.savefig(figname + '.png', dpi=600)
    # print 'Saved'
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
    st1.resample(50.0)
    st2.resample(50.0)
    print(st1)
    print(st2)
    ylist, x, comp_list = IntervalStackXCorr(st1,st2)
    print ylist
    print x
    print comp_list
    st2.resample(40.0)
    ylist2, x2, comp_list2 = IntervalStackXCorr(st1,st2)
    print ylist2
    print x2
    print comp_list2
    # confirm the error is zero
    maxdiff = np.max(np.abs(ylist2[0]-ylist[0]))
    print "Max diff = " + str(maxdiff)
    assert maxdiff < 0.05, "The xcor after resampling is not the same"
    saveXCorrPlot(ylist2,x2,'/g/data/ha3/','out2xcor',comp_list2)
    
    
