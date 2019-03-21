#!/usr/bin/env python
# coding: utf-8

import os
import sys

import numpy as np
import scipy
import matplotlib.dates
import matplotlib.pyplot as plt
import datetime
import dateutil
from dateutil import rrule

from seismic.ASDFdatabase import FederatedASDFDataSet

import obspy
from analytic_plot_utils import distance
from netCDF4 import Dataset as NCDataset

# Imports for plotting
from textwrap import wrap
from scipy import signal
from matplotlib.backends.backend_pdf import PdfPages


# Get station codes from file name
def stationCodes(filename):
    path, fname = os.path.split(filename)
    parts = fname.split('.')
    sta1 = '.'.join(parts[0:2])
    sta2 = '.'.join(parts[2:4])
    return (sta1, sta2)


def stationCoords(federated_ds, code, datetime):
    ds = federated_ds
    net, sta = code.split('.')
    sta_records = ds.get_stations(datetime, obspy.UTCDateTime(datetime) + 3600*24*30, network=net, station=sta)
    z_records = [r for r in sta_records if r[3][1:3] == 'HZ']
    assert len(z_records) == 1, z_records
    z_record = z_records[0]
    return z_record[4:6]


def stationDistance(federated_ds, code1, code2, datetime):
    coords1 = stationCoords(federated_ds, code1, datetime)
    coords2 = stationCoords(federated_ds, code2, datetime)
    return distance(coords1, coords2)


def computeEstimatedClockCorrections(rcf, snr_mask, ccf_masked, x_lag, pcf_cutoff_threshold):
    # Make an initial estimate of the shift, and only mask out a row if the Pearson coefficient
    # is less than the threshold AFTER applying the shift. Otherwise we will be masking out
    # some of the most interesting regions where shifts occur.
    correction = []
    ccf_shifted = []
    for row in ccf_masked:
        if np.ma.is_masked(row) or rcf is None:
            correction.append(np.nan)
            ccf_shifted.append(np.array([np.nan]*ccf_masked.shape[1]))
            continue

        # The logic here needs to be expanded to allow for possible mirroring of the CCF
        # as well as a shift.
        c3 = scipy.signal.correlate(rcf, row, mode='same')
        c3 /= np.max(c3)
        peak_index = np.argmax(c3)
        shift_size = int(peak_index - len(c3)/2)
        row_shifted = np.roll(row, shift_size)
        # Zero the rolled in values
        if shift_size > 0:
            row_shifted[0:shift_size] = 0
        elif shift_size < 0:
            row_shifted[shift_size:] = 0
        ccf_shifted.append(row_shifted)
        # Store first order corrections.
        peak_lag = x_lag[peak_index]
        correction.append(peak_lag)
    # end for
    correction = np.array(correction)
    ccf_shifted = np.array(ccf_shifted)
    # Recompute the RCF with first order clock corrections.
    rcf_corrected = np.nanmean(ccf_shifted[snr_mask, :], axis=0)

    # For the Pearson coeff threshold, apply it against the CORRECTED RCF after the application
    # of estimated clock corrections.
    row_rcf_crosscorr = []
    for i, row in enumerate(ccf_masked):
        if np.ma.is_masked(row) or rcf is None:
            row_rcf_crosscorr.append([np.nan]*ccf_masked.shape[1])
            continue

        pcf_corrected, _ = scipy.stats.pearsonr(rcf_corrected, ccf_shifted[i, :])
        if pcf_corrected < pcf_cutoff_threshold:
            correction[i] = np.nan
            row_rcf_crosscorr.append([np.nan]*ccf_masked.shape[1])
            continue
        # Compute second order corrections based on first order corrected RCF
        c3 = scipy.signal.correlate(rcf_corrected, row, mode='same')
        c3 /= np.max(c3)
        peak_index = np.argmax(c3)
        peak_lag = x_lag[peak_index]
        correction[i] = peak_lag
        row_rcf_crosscorr.append(c3)
    #end for
    row_rcf_crosscorr = np.array(row_rcf_crosscorr)

    return rcf_corrected, correction, row_rcf_crosscorr


def plotXcorrTimeseries(ax, x_lag, y_times, xcorr_data, use_formatter=False):

    np_times = np.array([datetime.datetime.utcfromtimestamp(v) for v in y_times]).astype('datetime64[s]')
    gx, gy = np.meshgrid(x_lag, np_times)
    im = ax.pcolormesh(gx, gy, xcorr_data, cmap='RdYlBu_r', vmin=0, vmax=1, rasterized=True)

    if use_formatter:
        date_formatter = matplotlib.dates.DateFormatter("%Y-%m-%d")
        date_locator = matplotlib.dates.WeekdayLocator(byweekday=rrule.SU)
        ax.yaxis.set_major_formatter(date_formatter)
        ax.yaxis.set_major_locator(date_locator)
    else:
        labels = np.datetime_as_string(np_times, unit='D')
        ax.set_yticks(np_times[::7])
        ax.set_yticklabels(labels[::7])

    ax.set_xlabel('Lag [s]')
    ax.set_ylabel('Days')
    
    ax_pos = ax.get_position()
    cax = plt.axes([ax_pos.x0 + 0.025, ax_pos.y1 - 0.1, 0.015, 0.08])

    plt.colorbar(im, cax=cax, orientation='vertical', ticks=[0, 1])


def plotRCF(ax, x_lag, rcf, rcf_corrected, snr_threshold):
    if rcf is not None:
        ax.axvline(x_lag[np.argmax(rcf)], c='#c66da9', lw=2,
                    label = '{:5.2f} s'.format(x_lag[np.argmax(rcf)]))
        ax.plot(x_lag, rcf, c='#42b3f4', 
                label="Reference CCF (RCF)\n"
                      "Based on Subset\n"
                      "with SNR > {}".format(snr_threshold))
        ax.plot(x_lag, rcf_corrected, '--', c='#00b75b', alpha=0.8, label="First order\ncorrected RCF")
        ax.legend()
    else:
        ax.text(0.5, 0.5, 'REFERENCE CCF:\nINSUFFICIENT SNR', horizontalalignment='center', 
         verticalalignment='center', transform=ax.transAxes, fontsize=16)

    ax.set_xticklabels([])


def plotStackedWindowCount(ax, x_nsw, y_times):
    ax.plot(x_nsw, y_times, c='#529664')
    ax.set_ylim((min(y_times), max(y_times)))
    ax.set_yticklabels([])
    ax.set_xlabel('\n'.join(wrap('# of Hourly Stacked Windows', 12)))
    xtl = ax.get_xticklabels()
    xtl[0].set_visible(False)
    xtl[-1].set_visible(False)


def plotSNRHistogram(ax, snr, time_window, nbins=10):
    ax.hist(snr.compressed(), fc='#42b3f4', ec='none', bins=nbins)
    ax.set_xlabel('SNR: Daily CCFs [-%d, %d]s'%(time_window, time_window))
    ax.set_ylabel('Frequency')
    xtl = ax.get_xticklabels()
    xtl[0].set_visible(False)
    xtl[-1].set_visible(False)


def plotPearsonCorrCoefficient(ax, rcf, ccf_masked, y_times):
    cc = []
    for row in ccf_masked:
        if np.ma.is_masked(row) or rcf is None:
            cc.append(np.nan)
        else:
            pcf, _ = scipy.stats.pearsonr(rcf, row)
            cc.append(pcf)
    # end for
    pcf = np.array(cc)
    # Compute CC mean
    ccav = np.mean(np.ma.masked_array(pcf, mask=np.isnan(pcf)))

    ax.plot(pcf, y_times, c='#d37f26')
    ax.set_ylim((min(y_times), max(y_times)))
    ax.set_yticklabels([])
    ax.set_xticks([0,1])
    ax.set_xlabel('\n'.join(wrap('Raw Pearson\nCoeff. (vs RCF)', 15)))
    ax.text(0.5, 0.98, '$PC_{av}$' + '={:3.3f}'.format(ccav), horizontalalignment='center', \
            verticalalignment='top', transform=ax.transAxes)
    ax.axvline(ccav, c='#d37f26', linestyle='--', lw=2, linewidth=1, alpha=0.7)


def plotEstimatedTimeshift(ax, x_lag, y_times, correction, annotation=None, row_rcf_crosscorr=None):

    if row_rcf_crosscorr is not None:
        # Line plot laid over the top of RCF * CCF
        np_times = np.array([datetime.datetime.utcfromtimestamp(v) for v in y_times]).astype('datetime64[s]')
        gx, gy = np.meshgrid(x_lag, np_times)
        plot_data = row_rcf_crosscorr
        crange_floor = 0.7
        plot_data[np.isnan(plot_data)] = crange_floor
        plot_data[(plot_data < crange_floor)] = crange_floor
        ax.pcolormesh(gx, gy, plot_data, cmap='RdYlBu_r', rasterized=True)
        ax.set_ylim((min(np_times), max(np_times)))
        xrange = 1.2*np.nanmax(np.abs(correction))
        ax.set_xlim((-xrange, xrange))
        col = '#ffffff'
        p = ax.plot(correction, np_times, 'o--', color=col, fillstyle='none', markersize=4, alpha=0.8)
#         ax.legend(p, ['Est. clock\ncorrection'], loc=1)
    else:
        # Plain line plot
        ax.plot(correction, y_times, 'o-', c='#f22e62', lw=1.5, fillstyle='none', markersize=4)
        ax.set_ylim((min(y_times), max(y_times)))
        ax.grid(":", color='#80808080')

    ax.set_yticklabels([])
    xtl = ax.get_xticklabels()
    xtl[0].set_visible(False)
    xtl[-1].set_visible(False)
    ax.set_xlabel('\n'.join(wrap('Estimated Timeshift [s]: RCF * CCF', 15)))
    if annotation is not None:
        ax.text(0.05, 0.98, annotation, color='#000000', horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes, rotation=90)


def plotXcorrFile(src_file, asdf_dataset, time_window, snr_threshold, show=True, underlay_rcf_xcorr=False,
                  pdf_file=None, png_file=None):
    # Read xcorr data
    xcdata = NCDataset(src_file, 'r')

    xc_start_times = xcdata.variables['IntervalStartTimes'][:] # sTimes
    xc_end_times = xcdata.variables['IntervalEndTimes'][:] # eTimes
    xc_lag = xcdata.variables['lag'][:] # lag
    xc_xcorr = xcdata.variables['xcorr'][:, :] # xcorr
    xc_nStackedWindows = xcdata.variables['NumStackedWindows'][:] # nStackedWindows
    xcdata.close()
    xcdata = None

    start_utc_time = obspy.UTCDateTime(xc_start_times[0])
    end_utc_time = obspy.UTCDateTime(xc_end_times[-1])

    start_time = str(start_utc_time)
    end_time = str(end_utc_time)

    origin_code, dest_code = stationCodes(src_file)
    dist = stationDistance(asdf_dataset, origin_code, dest_code, start_time)
    
    # Extract primary data
    lagIndices = np.squeeze(np.argwhere(np.fabs(np.round(xc_lag, decimals=2)) == time_window))
    sTimes = xc_start_times
    lag = xc_lag[lagIndices[0]:lagIndices[1]]
    ccf = xc_xcorr[:, lagIndices[0]:lagIndices[1]]
    nsw = xc_nStackedWindows

    # Compute derived quantities used by multiple axes
    zero_row_mask = (np.all(ccf == 0, axis=1))
    valid_mask = np.ones_like(ccf)
    valid_mask[zero_row_mask, :] = 0
    valid_mask = (valid_mask > 0)
    ccfMasked = np.ma.masked_array(ccf, mask=~valid_mask)
    snr = np.nanmax(ccfMasked, axis=1) / np.nanstd(ccfMasked, axis=1)
    if np.any(snr > snr_threshold):
        snr_mask = (snr > snr_threshold)
        rcf = np.nanmean(ccfMasked[snr_mask, :], axis=0)
    else:
        snr_mask = None
        rcf = None

    PCF_CUTOFF_THRESHOLD = 0.5
    rcf_corrected, correction, row_rcf_crosscorr = \
        computeEstimatedClockCorrections(rcf, snr_mask, ccfMasked, lag, PCF_CUTOFF_THRESHOLD)

    #-----------------------------------------------------------------------
    # Master layout and plotting code
    
    fig = plt.figure(figsize=(11.69,16.53))
    fig.suptitle("Station: {}, Dist. to {}: {:3.2f} km".format(origin_code, dest_code, dist), fontsize = 16, y=1)

    ax1 = fig.add_axes([0.1, 0.075, 0.5, 0.725])

    labelPad = 0.05
    ax2 = fig.add_axes([0.1, 0.8, 0.5, 0.175]) # reference CCF (accumulation of daily CCFs)
    ax3 = fig.add_axes([0.6, 0.075, 0.1, 0.725]) # number of stacked windows
    ax4 = fig.add_axes([0.6 + labelPad, 0.8 + 0.6*labelPad, 0.345, 0.175 - 0.6*labelPad]) # SNR histogram
    ax5 = fig.add_axes([0.7, 0.075, 0.1, 0.725]) # Pearson coeff
    ax6 = fig.add_axes([0.8, 0.075, 0.195, 0.725]) # estimate timeshifts

    # Plot CCF image =======================
    plotXcorrTimeseries(ax1, lag, sTimes, ccf)

    # Plot CCF-template (reference CCF) ===========
    plotRCF(ax2, lag, rcf, rcf_corrected, snr_threshold)

    # Plot number of stacked windows ==============
    plotStackedWindowCount(ax3, nsw, sTimes)

    # Plot histogram
    plotSNRHistogram(ax4, snr, time_window)

    # Plot Pearson correlation coefficient=========
    plotPearsonCorrCoefficient(ax5, rcf, ccfMasked, sTimes)

    # plot Timeshift =====================
    annotation = 'Min. corrected Pearson Coeff={:3.3f}'.format(PCF_CUTOFF_THRESHOLD)
    if underlay_rcf_xcorr:
        plotEstimatedTimeshift(ax6, lag, sTimes, correction, annotation=annotation, 
                               row_rcf_crosscorr=row_rcf_crosscorr)
    else:
        plotEstimatedTimeshift(ax6, lag, sTimes, correction, annotation=annotation)

    # Print and display
    if pdf_file is not None:
        pdf_out = PdfPages(pdf_file)
        pdf_out.savefig(plt.gcf(), dpi=600)
        pdf_out.close()

    if png_file is not None:
        plt.savefig(png_file, dpi=150)

    if show:
        plt.show()

    plt.close()

