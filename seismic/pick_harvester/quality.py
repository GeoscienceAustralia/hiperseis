"""
Description:
    Computes quality measures of picks using various methods
References:

CreationDate:   24/01/19
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     24/01/19   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import glob, os, sys
import pywt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from numpy import unravel_index
import numpy as np
import traceback

def compute_quality_measures(trc, scales, plotinfo=None):
    # function for a least-squares line fit
    def function(x, A, B):
        return A * x + B
    # end func

    try:
        # =======================================
        # Compute wavelets-based quality estimate
        # =======================================
        cwt, freqs = pywt.cwt(trc, scales, 'gaus8', trc.stats.delta)
        ps = np.fabs(cwt) ** 2
        ps = pywt.threshold(ps, np.std(ps), mode='soft', substitute=1)

        psbefore   = ps[:, :ps.shape[1] / 2]
        psafter    = ps[:, ps.shape[1] / 2:]
        idx_max_before = unravel_index(psbefore.argmax(), psbefore.shape)
        idx_max_after  = unravel_index(psafter.argmax(), psafter.shape)

        #before = psbefore[idx_max_after[0], :]
        before = np.amax(psbefore, axis=0)
        after  = np.amax(psafter, axis=0)
        cwtsnr = np.mean(after[after > np.std(after)]) / \
                 np.mean(before)

        # =======================================
        # Compute slope-based quality estimate
        # =======================================
        times = trc.times() - trc.times().max() / 2

        mid = len(trc.times()) / 2
        timesa = times[mid:]
        timesb = times[:mid]
        ab = np.cumsum(np.fabs(trc.data))

        popta, _ = curve_fit(function, timesa, ab[mid:])
        poptb, _ = curve_fit(function, timesb, ab[:mid])
        slope_ratio = popta[0] / poptb[0]

        # =======================================
        # Plot and save results
        # =======================================
        if(plotinfo):
            ei = plotinfo['eventid']
            ot = str(plotinfo['origintime'])
            mag = plotinfo['mag']
            net = plotinfo['net']
            sta = plotinfo['sta']
            phase = plotinfo['phase']
            ppsnr = plotinfo['ppsnr']
            pickid = plotinfo['pickid']
            of = plotinfo['outputfolder']

            fig, axes = plt.subplots(ncols=1, nrows=3)

            fig.suptitle("""OriginTime: %s\nMagnitude: %4.2f
            PhasePapySNR: %4.2f
            WaveletsSNR: %4.2f
            SlopeRatio: %4.2f""" % (ot, mag, ppsnr, cwtsnr, slope_ratio), fontsize=6)

            x = times
            y = freqs
            x, y = np.meshgrid(x, y)

            axes[0].plot(times, trc.data)
            axes[1].pcolormesh(x, y, ps, cmap='jet')
            axes[2].plot(timesa, ab[mid:], c='r')
            axes[2].plot(timesb, ab[:mid], c='b')
            axes[2].plot(timesa, function(timesa, *popta), c='k', linewidth=0.5, linestyle='--')
            axes[2].plot(timesb, function(timesb, *poptb), c='k', linewidth=0.5, linestyle='--')

            axes[0].set_xlim(times.min(), times.max())
            axes[0].tick_params(labelbottom=False)
            axes[1].set_ylabel('Freq [Hz]')
            axes[1].tick_params(labelbottom=False)
            axes[2].set_xlim(times.min(), times.max())
            axes[2].set_xlabel('Time [s]')
            axes[2].set_ylabel('Cumm. Length')

            ofn = os.path.join(of, '%s.%d.%s.%s.%d.png' % (phase, ei, net, sta, pickid))
            plt.savefig(ofn, dpi=300)
            plt.close()
        # end if
    except:
        traceback.print_exc()
    # end try

    return cwtsnr, slope_ratio
# end func