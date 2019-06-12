#!/usr/bin/env python
# coding: utf-8
"""
Functions for computing estimated GPS clock corrections based on station pair cross-correlation
and plotting in standard layout.
"""

# pylint: disable=invalid-name, no-name-in-module, broad-except

import os
import sys
import datetime
import time
import math
import gc
import glob

from textwrap import wrap

import numpy as np
import scipy
from scipy import signal
import pandas as pd
import matplotlib.dates
import matplotlib.pyplot as plt
from dateutil import rrule
import click
from sklearn.cluster import dbscan
from scipy.interpolate import LSQUnivariateSpline
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

# import obspy
from netCDF4 import Dataset as NCDataset
from tqdm.auto import tqdm

from seismic.ASDFdatabase import FederatedASDFDataSet
from seismic.xcorqc.analytic_plot_utils import distance, timestamps_to_plottable_datetimes


class XcorrClockAnalyzer:
    """
    Helper class for bundling of preprocessed cross-correlation data before plotting
    or subsequent processing.
    """
    def __init__(self, src_file, time_window, snr_threshold, pcf_cutoff_threshold):
        """Constructor

        :param src_file: Source .nc file
        :type src_file: str
        :param time_window: Time window to analyze (maximum detectable time shift), in seconds
        :type time_window: int
        :param snr_threshold: Minimum SNR threshold for a given day to include in stacked correlation function (RCF)
        :type snr_threshold: int
        :param pcf_cutoff_threshold: Minimum Pearson correlation factor (RCF * CCF) to accept for generating a time
            shift estimate for a given day, between 0.0 and 1.0.
        :type pcf_cutoff_threshold: float
        """
        # INPUTS
        # Input file
        self.src_file = src_file
        # Size of lag window to analyse
        self.time_window = time_window
        # Minimum SNR per sample
        self.snr_threshold = snr_threshold
        # Minimum Pearson correlation factor to accept for a given day's CCF
        assert pcf_cutoff_threshold >= 0.0 and pcf_cutoff_threshold <= 1.0
        self.pcf_cutoff_threshold = pcf_cutoff_threshold

        # OUTPUTS
        # Reference correlation function (RCF) as per Hable et. al (2018)
        # :type rcf: np.array
        self.rcf = None
        # Timestamps corresponding to the time series data in ccf, lag, nsw
        self.start_times = None
        # Cross-correlation function sample
        self.ccf = None
        # Masked ccf according to quality criteria
        # :type ccf_masked: numpy.ma.masked_array
        self.ccf_masked = None
        # Signal to noise ratio
        self.snr = None
        # Mask of which samples meet minimum SNR criteria
        # :type snr_mask: numpy array mask
        self.snr_mask = None
        # Horizontal axis values (lag window)
        # :type lag: numpy.array
        self.lag = None
        # Number of stacked windows per sample
        self.nsw = None

        # Generate outputs
        self._preprocess()
        self._compute_estimated_clock_corrections()

        # Get analytical time series usable for clustering
        nan_mask = np.isnan(self.raw_correction)
        self.correction_times_clean = self.float_start_times[~nan_mask]
        self.corrections_clean = self.raw_correction[~nan_mask]
        # Compute a median filtered slope metric
        grad = np.gradient(self.corrections_clean, self.correction_times_clean, edge_order=1)
        grad_med5 = signal.medfilt(grad, 5)
        self.corrections_slope = grad_med5
    # end func

    def get_corrections_time_series(self):
        return self.correction_times_clean, self.corrections_clean

    def do_clustering(self, coeffs):
        """
        Do DBSCAN clustering on the corrections.

        :param coeffs: Triplet of distance coefficients, corresponding to the sensitivity of the clustering to point
            separation along 1) x-axis (time), 2) y-axis (correction) and 3) slope (drift rate)
        :type coeffs: tuple(float, float, float)
        :return: Results of sklearn.cluster.dbscan (refer to third party documentation)
        """
        sec_per_week = 7 * 24 * 3600
        def _temporalDist2DSlope(p0, p1, coeffs):
            return math.sqrt((coeffs[0] * (p1[0] - p0[0])) ** 2 +
                             (coeffs[1] * sec_per_week * (p1[1] - p0[1])) ** 2 +
                             (coeffs[2] * sec_per_week * sec_per_week * (p1[2] - p0[2])) ** 2)

        data = np.column_stack((self.correction_times_clean, self.corrections_clean, self.corrections_slope))
        ind, ids = dbscan(data, eps=2 * sec_per_week, min_samples=7,
                          metric=lambda p0, p1: _temporalDist2DSlope(p0, p1, coeffs))
        return ind, ids
    # end func

    def do_spline_regression(self, group_ids, regression_degree):
        """
        Do univariate spline regression on each cluster of points.

        :param group_ids: Cluster IDS generated from do_clustering()
        :param regression_degree: Desired degree of curve fit for each cluster, one for each non-negative cluster ID
        :return: dict of regressors that can be applied to arbitrary time values
        """
        assert isinstance(regression_degree, dict)
        regressions = {}
        gids = set(group_ids[group_ids != -1])
        for i in gids:
            degree = regression_degree[i]
            mask_group = (group_ids == i)
            # Perform regression
            x = self.correction_times_clean[mask_group]
            y = self.corrections_clean[mask_group]
            # Constrain the knot frequency to be proportional to the order, so that we don't overfit
            # to high frequency variations in the data just for the sake of lowering the spline residual.
            t = np.linspace(x[0], x[-1], degree - 1)
            r = LSQUnivariateSpline(x, y, t[1:-1], k=degree)
            regressions[i] = r

        return regressions
    # end func

    def do_spline_resampling(self, group_ids, regressors, sampling_period_seconds):
        """Using pre-computed regressors, resample every cluster at a prescribed frequency

        :param group_ids:
        :param regressors:
        :param sampling_period_seconds:
        :return:
        """
        # Dict of daily spaced time values and computed correction, since source data time
        # points might not be uniformly distributed. Keyed by group ID.
        regular_corrections = {}
        gids = set(group_ids[group_ids != -1])
        for i in gids:
            mask_group = (group_ids == i)
            # Generate uniform daily times at which to compute corrections
            x = self.correction_times_clean[mask_group]
            timestamp_min = min(x)
            timestamp_max = max(x)
            num_samples = np.round((timestamp_max - timestamp_min) / sampling_period_seconds)
            lin_times = np.linspace(timestamp_min, timestamp_max, num_samples + 1)
            lin_corrections = regressors[i](lin_times)
            regular_corrections[i] = {'times': lin_times, 'corrections': lin_corrections}

        return regular_corrections
    # end func

    def _preprocess(self):
        xcdata = NCDataset(self.src_file, 'r')

        xc_start_times = xcdata.variables['IntervalStartTimes'][:]
        # xc_end_times = xcdata.variables['IntervalEndTimes'][:]
        xc_lag = xcdata.variables['lag'][:]
        xc_xcorr = xcdata.variables['xcorr'][:, :]
        xc_num_stacked_windows = xcdata.variables['NumStackedWindows'][:]
        xcdata.close()
        xcdata = None

        # start_utc_time = obspy.UTCDateTime(xc_start_times[0])
        # end_utc_time = obspy.UTCDateTime(xc_end_times[-1])
        # start_time = str(start_utc_time)
        # end_time = str(end_utc_time)
        # print("Date range for file {}:\n    {} -- {}".format(self.src_file, start_time, end_time))

        # Extract primary data
        lag_indices = np.squeeze(np.argwhere(np.fabs(np.round(xc_lag, decimals=2)) <= self.time_window))
        self.start_times = xc_start_times
        self.lag = xc_lag[lag_indices[0]:lag_indices[-1]]
        self.ccf = xc_xcorr[:, lag_indices[0]:lag_indices[-1]]
        self.nsw = xc_num_stacked_windows

        # Compute derived quantities used by multiple axes
        zero_row_mask = (np.all(self.ccf == 0, axis=1))
        valid_mask = np.ones_like(self.ccf)
        valid_mask[zero_row_mask, :] = 0
        valid_mask = (valid_mask > 0)
        self.ccf_masked = np.ma.masked_array(self.ccf, mask=~valid_mask)
        self.snr = np.nanmax(self.ccf_masked, axis=1) / np.nanstd(self.ccf_masked, axis=1)
        if np.any(self.snr > self.snr_threshold):
            self.snr_mask = (self.snr > self.snr_threshold)
            self.rcf = np.nanmean(self.ccf_masked[self.snr_mask, :], axis=0)
        self.float_start_times = np.array([float(v) for v in self.start_times])

    # end func

    def _compute_estimated_clock_corrections(self):
        """
        Compute the estimated GPS clock corrections given series of cross-correlation functions
        and an overall reference correlation function (rcf, the mean of valid samples of the
        cross-correlation time series).

        :return: Corrected RCF after first order corrections; raw estimated clock corrections;
                 final per-sample cross-correlation between corrected RCF and each sample
        :rtype: numpy.array [1D]
        """
        # Make an initial estimate of the shift, and only mask out a row if the Pearson coefficient
        # is less than the threshold AFTER applying the shift. Otherwise we will be masking out
        # some of the most interesting regions where shifts occur.
        raw_correction = []
        ccf_shifted = []
        for row in self.ccf_masked:
            if np.ma.is_masked(row) or self.rcf is None:
                raw_correction.append(np.nan)
                ccf_shifted.append(np.array([np.nan] * self.ccf_masked.shape[1]))
                continue

            # The logic here needs to be expanded to allow for possible mirroring of the CCF
            # as well as a shift.
            c3 = scipy.signal.correlate(self.rcf, row, mode='same')
            c3 /= np.max(c3)
            peak_index = np.argmax(c3)
            shift_size = int(peak_index - len(c3) / 2)
            row_shifted = np.roll(row, shift_size)
            # Zero the rolled in values
            if shift_size > 0:
                row_shifted[0:shift_size] = 0
            elif shift_size < 0:
                row_shifted[shift_size:] = 0
            ccf_shifted.append(row_shifted)
            # Store first order corrections.
            peak_lag = self.lag[peak_index]
            raw_correction.append(peak_lag)
        # end for
        raw_correction = np.array(raw_correction)
        ccf_shifted = np.array(ccf_shifted)
        # Recompute the RCF with first order clock corrections.
        rcf_corrected = np.nanmean(ccf_shifted[self.snr_mask, :], axis=0)

        # For the Pearson coeff threshold, apply it against the CORRECTED RCF after the application
        # of estimated clock corrections.
        row_rcf_crosscorr = []
        for i, row in enumerate(self.ccf_masked):
            if np.ma.is_masked(row) or self.rcf is None:
                row_rcf_crosscorr.append([np.nan] * self.ccf_masked.shape[1])
                continue

            pcf_corrected, _ = scipy.stats.pearsonr(rcf_corrected, ccf_shifted[i, :])
            if pcf_corrected < self.pcf_cutoff_threshold:
                raw_correction[i] = np.nan
                row_rcf_crosscorr.append([np.nan] * self.ccf_masked.shape[1])
                continue
            # Compute second order corrections based on first order corrected RCF
            c3 = scipy.signal.correlate(rcf_corrected, row, mode='same')
            c3 /= np.max(c3)
            peak_index = np.argmax(c3)
            peak_lag = self.lag[peak_index]
            raw_correction[i] = peak_lag
            row_rcf_crosscorr.append(c3)
        # end for
        row_rcf_crosscorr = np.array(row_rcf_crosscorr)

        self.raw_correction = raw_correction
        self.rcf_corrected = rcf_corrected
        self.row_rcf_crosscorr = row_rcf_crosscorr
    # end func

    def plot_clusters(self, ax, ids, coeffs, stn_code=''):
        """Plot the distinct clusters color coded by cluster ID, with underlying corrections shown in gray.

        :param ax:
        :param ids:
        :param coeffs:
        :return:
        """
        plot_times = timestamps_to_plottable_datetimes(self.correction_times_clean)
        ax.plot(plot_times, self.corrections_clean, 'x', color="#808080")
        for i in range(max(ids) + 1):
            mask_group = (ids == i)
            ax.plot(plot_times[mask_group], self.corrections_clean[mask_group], 'o-',
                    color='C{}'.format(i), markersize=6, fillstyle='none', alpha=0.7)

        time_formatter = matplotlib.dates.DateFormatter("%Y-%m-%d")
        ax.xaxis.set_major_formatter(time_formatter)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.grid(':', color="#808080", zorder=0, alpha=0.5)
        ax.set_xlabel('Day', fontsize=14)
        ax.set_ylabel('Correction (sec)', fontsize=14)
        ax.text(0.02, 0.96, "Cluster coeffs: {}".format(coeffs), fontsize=12, color="#606060",
                transform=ax.transAxes)
        if stn_code:
            ax.set_title("Station {} first order corrections groups".format(stn_code), fontsize=20)
    # end func

    def plot_regressors(self, ax, ids, regressors, stn_code=''):
        """Plot regressor functions on top of original data

        :param ax:
        :param ids:
        :param regressors:
        :param stn_code:
        :return:
        """
        assert isinstance(regressors, dict)
        gids = set(ids[ids != -1])
        # Corrections from regression function fit
        correction_fit = np.zeros_like(self.corrections_clean)
        correction_fit[:] = np.nan
        for i in gids:
            mask_group = (ids == i)
            # Apply regression
            x = self.correction_times_clean[mask_group]
            # Compute fitted values
            correction_fit[mask_group] = regressors[i](x)

        plot_times = timestamps_to_plottable_datetimes(self.correction_times_clean)
        for i in gids:
            mask_group = (ids == i)
            ax.plot(plot_times[mask_group], self.corrections_clean[mask_group], 'o', color='C{}'.format(i),
                     markersize=5, fillstyle='none')
            ax.plot(plot_times[mask_group], correction_fit[mask_group], color='C{}'.format(i))

        time_formatter = matplotlib.dates.DateFormatter("%Y-%m-%d")
        ax.xaxis.set_major_formatter(time_formatter)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.grid(':', color="#808080", zorder=0, alpha=0.5)
        ax.set_xlabel('Day', fontsize=14)
        ax.set_ylabel('Correction (sec)', fontsize=14)
        if stn_code:
            ax.set_title("Station {} corrections groups with regressions to sample times".format(stn_code), fontsize=20)
    # end func

    def plot_resampled_clusters(self, ax, ids, resampled_corrections, stn_code=''): # ax, ids, regular_corrections, FULL_CODE
        """Plot resampled regressor functions per cluster on top of original data

        :param ax:
        :param ids:
        :param regressors:
        :param stn_code:
        :return:
        """
        assert isinstance(resampled_corrections, dict)
        gids = set(ids[ids != -1])
        plot_times = timestamps_to_plottable_datetimes(self.correction_times_clean)
        for i in gids:
            mask_group = (ids == i)
            ax.plot(plot_times[mask_group], self.corrections_clean[mask_group], 'o', color='#808080'.format(i),
                    markersize=5, fillstyle='none')
            plt_times_i = timestamps_to_plottable_datetimes(resampled_corrections[i]['times'])
            ax.plot(plt_times_i, resampled_corrections[i]['corrections'], 'o', color='C{}'.format(i), fillstyle='none')

        time_formatter = matplotlib.dates.DateFormatter("%Y-%m-%d")
        ax.xaxis.set_major_formatter(time_formatter)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.grid(':', color="#808080", zorder=0, alpha=0.5)
        ax.set_xlabel('Day', fontsize=14)
        ax.set_ylabel('Correction (sec)', fontsize=14)
        if stn_code:
            ax.set_title("Station {} corrections groups with regressions to daily samples".format(stn_code), fontsize=20)
    # end func

#end class

def station_codes(filename):
    """
    Convert a netCDF4 .nc filename generated by correlator to the corresponding
    station codes in the format ``NETWORK.STATION``

    Assumed format: ``NET1.STAT1.NET2.STA2.*.nc``

    :param filename: The ``.nc`` filename from which to extract the station and network codes
    :type filename: str
    :return: Station code for each station (code1, code2)
    :rtype: tuple(str, str)
    """
    _, fname = os.path.split(filename)
    parts = fname.split('.')
    sta1 = '.'.join(parts[0:2])
    sta2 = '.'.join(parts[2:4])
    return (sta1, sta2)


def station_distance(federated_ds, code1, code2):
    """
    Determine the distance in km between a pair of stations at a given time.

    :param federated_ds: Federated dataset to query for station coordinates
    :type federated_ds: seismic.ASDFdatabase.FederatedASDFDataSet
    :param code1: Station and network code in the format ``NETWORK.STATION`` for first station
    :type code1: str
    :param code2: Station and network code in the format ``NETWORK.STATION`` for second station
    :type code2: str
    :return: Distance between stations in kilometers
    :rtype: float
    """
    coords1 = federated_ds.unique_coordinates[code1]
    coords2 = federated_ds.unique_coordinates[code2]
    return distance(coords1, coords2)


def plot_xcorr_time_series(ax, x_lag, y_times, xcorr_data, use_formatter=False):

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


def plot_reference_correlation_function(ax, x_lag, rcf, rcf_corrected, snr_threshold):
    if rcf is not None:
        ax.axvline(x_lag[np.argmax(rcf)], c='#c66da9', lw=2,
                   label='{:5.2f} s'.format(x_lag[np.argmax(rcf)]))
        ax.plot(x_lag, rcf, c='#42b3f4',
                label="Reference CCF (RCF)\n"
                      "Based on Subset\n"
                      "with SNR > {}".format(snr_threshold))
        ax.plot(x_lag, rcf_corrected, '--', c='#00b75b', alpha=0.8, label="First order\ncorrected RCF")
        ax.legend(loc=1)
    else:
        ax.text(0.5, 0.5, 'REFERENCE CCF:\nINSUFFICIENT SNR', horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes, fontsize=16)

    ax.set_xticklabels([])


def plot_stacked_window_count(ax, x_nsw, y_times):
    ax.plot(x_nsw, y_times, c='#529664')
    ax.set_ylim((min(y_times), max(y_times)))
    ax.set_yticklabels([])
    ax.set_xlabel('\n'.join(wrap('# of Hourly Stacked Windows', 12)))
    xtl = ax.get_xticklabels()
    xtl[0].set_visible(False)
    xtl[-1].set_visible(False)


def plot_snr_histogram(ax, snr, time_window, nbins=10):
    ax.hist(snr.compressed(), fc='#42b3f4', ec='none', bins=nbins)
    ax.set_xlabel('SNR: Daily CCFs [-{0:d}, {0:d}]s'.format(time_window))
    ax.set_ylabel('Frequency')
    xtl = ax.get_xticklabels()
    xtl[0].set_visible(False)
    xtl[-1].set_visible(False)


def plot_pearson_corr_coeff(ax, rcf, ccf_masked, y_times):
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
    ax.set_xticks([0, 1])
    ax.set_xlabel('\n'.join(wrap('Raw Pearson\nCoeff. (vs RCF)', 15)))
    ax.text(0.5, 0.98, '$PC_{av}$' + '={:3.3f}'.format(ccav), horizontalalignment='center',
            verticalalignment='top', transform=ax.transAxes)
    ax.axvline(ccav, c='#d37f26', linestyle='--', lw=2, linewidth=1, alpha=0.7)


def plot_estimated_timeshift(ax, x_lag, y_times, correction, annotation=None, row_rcf_crosscorr=None):

    if row_rcf_crosscorr is not None:
        # Line plot laid over the top of RCF * CCF
        np_times = np.array([datetime.datetime.utcfromtimestamp(v) for v in y_times]).astype('datetime64[s]')
        gx, gy = np.meshgrid(x_lag, np_times)
        plot_data = row_rcf_crosscorr
        crange_floor = 0.7
        plot_data[np.isnan(plot_data)] = crange_floor
        plot_data[(plot_data < crange_floor)] = crange_floor
        nan_mask = np.isnan(correction)
        plot_data[nan_mask, :] = np.nan
        ax.pcolormesh(gx, gy, plot_data, cmap='RdYlBu_r', rasterized=True)
        ax.set_ylim((min(np_times), max(np_times)))
        xrange = 1.2 * np.nanmax(np.abs(correction))
        ax.set_xlim((-xrange, xrange))
        col = '#ffffff'
        ax.plot(correction, np_times, 'o--', color=col, fillstyle='none', markersize=4, alpha=0.8)
    else:
        # Plain line plot
        np_times = np.array([datetime.datetime.utcfromtimestamp(v) for v in y_times]).astype('datetime64[s]')
        ax.plot(correction, np_times, 'o-', c='#f22e62', lw=1.5, fillstyle='none', markersize=4)
        ax.set_ylim((min(np_times), max(np_times)))
        xlim = list(ax.get_xlim())
        xlim[0] = min(xlim[0], -1)
        xlim[1] = max(xlim[1], 1)
        ax.set_xlim(tuple(xlim))
        ax.grid(":", color='#80808080')

    ytl = ax.get_yticklabels()
    for y in ytl:
        y.set_visible(False)
    xtl = ax.get_xticklabels()
    xtl[0].set_visible(False)
    xtl[-1].set_visible(False)
    ax.set_xlabel('\n'.join(wrap('Estimated Timeshift [s]: RCF * CCF', 15)))
    if annotation is not None:
        ax.text(0.05, 0.98, annotation, color='#000000', horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes, rotation=90)


def plot_xcorr_file_clock_analysis(src_file, asdf_dataset, time_window, snr_threshold, pearson_correlation_factor,
                                   show=True, underlay_rcf_xcorr=False,
                                   pdf_file=None, png_file=None, title_tag='',
                                   settings=None):
    # Read and preprocess xcorr data
    xcorr_ca = XcorrClockAnalyzer(src_file, time_window, snr_threshold, pearson_correlation_factor)
    raw_correction = xcorr_ca.raw_correction
    rcf_corrected = xcorr_ca.rcf_corrected
    row_rcf_crosscorr = xcorr_ca.row_rcf_crosscorr

    origin_code, dest_code = station_codes(src_file)
    dist = station_distance(asdf_dataset, origin_code, dest_code)

    # -----------------------------------------------------------------------
    # Master layout and plotting code

    fig = plt.figure(figsize=(11.69, 16.53))
    fig.suptitle("Station: {}, Dist. to {}: {:3.2f} km{}".format(
        origin_code, dest_code, dist, '' if not title_tag else ' ' + title_tag), fontsize=16, y=1)

    ax1 = fig.add_axes([0.1, 0.075, 0.5, 0.725])

    label_pad = 0.05
    ax2 = fig.add_axes([0.1, 0.8, 0.5, 0.175])  # reference CCF (accumulation of daily CCFs)
    ax3 = fig.add_axes([0.6, 0.075, 0.1, 0.725])  # number of stacked windows
    ax4 = fig.add_axes([0.6 + label_pad, 0.8 + 0.6 * label_pad, 0.345, 0.175 - 0.6 * label_pad])  # SNR histogram
    ax5 = fig.add_axes([0.7, 0.075, 0.1, 0.725])  # Pearson coeff
    ax6 = fig.add_axes([0.8, 0.075, 0.195, 0.725])  # estimate timeshifts

    # Plot CCF image =======================
    plot_xcorr_time_series(ax1, xcorr_ca.lag, xcorr_ca.start_times, xcorr_ca.ccf)

    # Plot CCF-template (reference CCF) ===========
    plot_reference_correlation_function(ax2, xcorr_ca.lag, xcorr_ca.rcf, rcf_corrected, snr_threshold)

    # If settings are provided, plot them as text box overlay in ax2
    if settings is not None:
        def setstr(name):
            return '{}: {}'.format(name, settings[name])
        settings_str = '\n'.join([setstr('INTERVAL_SECONDS'), setstr('WINDOW_SECONDS'),
                                  setstr('--resample-rate'), setstr('--fmin'), setstr('--fmax'),
                                  setstr('--clip-to-2std'), setstr('--one-bit-normalize'),
                                  setstr('--envelope-normalize')])
        if '--whitening' in settings:
            settings_str += '\n' + setstr('--whitening')
        ax2.text(0.02, 0.97, settings_str, fontsize=8, alpha=0.8,
                 horizontalalignment='left', verticalalignment='top',
                 transform=ax2.transAxes, bbox=dict(fc='#e0e0e080', ec='#80808080'))

    # Plot number of stacked windows ==============
    plot_stacked_window_count(ax3, xcorr_ca.nsw, xcorr_ca.start_times)

    # Plot histogram
    plot_snr_histogram(ax4, xcorr_ca.snr, time_window)

    # Plot Pearson correlation coefficient=========
    plot_pearson_corr_coeff(ax5, xcorr_ca.rcf, xcorr_ca.ccf_masked, xcorr_ca.start_times)

    # plot Timeshift =====================
    annotation = 'Min. corrected Pearson Coeff={:3.3f}'.format(pearson_correlation_factor)
    if underlay_rcf_xcorr:
        plot_estimated_timeshift(ax6, xcorr_ca.lag, xcorr_ca.start_times, raw_correction, annotation=annotation,
                                 row_rcf_crosscorr=row_rcf_crosscorr)
    else:
        plot_estimated_timeshift(ax6, xcorr_ca.lag, xcorr_ca.start_times, raw_correction, annotation=annotation)

    # Print and display
    if pdf_file is not None:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf_out = PdfPages(pdf_file)
        pdf_out.savefig(plt.gcf(), dpi=600)
        pdf_out.close()

    if png_file is not None:
        plt.savefig(png_file, dpi=150)

    if show:
        plt.show()

    # If no output options, return the figure handle for further use by caller
    if pdf_file is None and png_file is None and not show:
        return xcorr_ca, fig
    else:
        # Need to clean all this up explicitly, otherwise huge memory leaks when run
        # from ipython notebook.
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        ax5.clear()
        ax6.clear()
        fig.clf()
        plt.close('all')
        return xcorr_ca, None


def read_correlator_config(nc_file):
    """
    Read the correlator settings used for given nc file.

    :param nc_file: File name of the .nc file containing the cross-correlation data.
    :type nc_file: str
    :return: Pandas Series with named fields whose values are the runtime settings used for the .nc file
    :rtype: pandas.Series
    """
    folder, fname = os.path.split(nc_file)
    base, _ = os.path.splitext(fname)
    base_parts = base.split('.')
    if len(base_parts) > 4:
        timestamp = '.'.join(base_parts[4:])
        config_filename = '.'.join(['correlator', timestamp, 'cfg'])
        config_file = os.path.join(folder, config_filename)
        if os.path.exists(config_file):
            settings_df = pd.read_csv(config_file, sep=':', names=['setting', 'value'],
                                      index_col=0, skiprows=1, skipinitialspace=True,
                                      squeeze=True, converters={'setting': str.strip})
            title_tag = '({}-{} Hz)'.format(settings_df['--fmin'], settings_df['--fmax'])
        else:
            print('WARNING: Settings file {} not found!'.format(config_file))
            settings_df = None
            title_tag = ''
    else:
        settings_df = None
        title_tag = ''

    return settings_df, title_tag


def batch_process_xcorr(src_files, dataset, time_window, snr_threshold, pearson_cutoff_factor=0.5, save_plots=True,
                        underlay_rcf_xcorr=False, force_save=False):
    """
    Process a batch of .nc files to generate standard visualization graphics. PNG files are output alongside the
    source .nc file. To suppress file output, set save_plots=False.

    :param src_files: List of files to process
    :type src_files: Iterable of str
    :param dataset: Dataset to be used to ascertain the distance between stations.
    :type dataset: FederatedASDFDataset
    :param time_window: Lag time window to plot (plus or minus this value in seconds)
    :type time_window: float
    :param snr_threshold: Minimum signal to noise ratio for samples to be included into the clock lag estimate
    :type snr_threshold: float
    :param save_plots: Whether to save plots to file, defaults to True
    :param save_plots: bool, optional
    :param underlay_rcf_xcorr: Show the individual correlation of row sample with RCF beneath the computed time lag,
        defaults to False
    :param underlay_rcf_xcorr: bool, optional
    :return: List of files for which processing failed, and associated error.
    :rtype: list(tuple(str, str))
    """
    PY2 = (sys.version_info[0] == 2)

    pbar = tqdm(total=len(src_files), dynamic_ncols=True)
    found_preexisting = False
    failed_files = []
    skipped_count = 0
    success_count = 0
    for src_file in src_files:
        _, base_file = os.path.split(src_file)
        pbar.set_description(base_file)
        # Sleep to ensure progress bar is refreshed
        time.sleep(0.2)

        if not os.path.exists(src_file):
            tqdm.write("ERROR! File {} not found!".format(src_file))
            failed_files.append((src_file, "File not found!"))
            continue

        # Extract timestamp from nc filename if available
        settings, title_tag = read_correlator_config(src_file)

        try:
            if save_plots:
                basename, _ = os.path.splitext(src_file)
                png_file = basename + ".png"
                # If png file already exists and has later timestamp than src_file, then skip it.
                if os.path.exists(png_file):
                    src_file_time = os.path.getmtime(src_file)
                    png_file_time = os.path.getmtime(png_file)
                    png_file_size = os.stat(png_file).st_size
                    if not force_save and (png_file_time > src_file_time) and (png_file_size > 0):
                        tqdm.write("PNG file {} is more recent than source file {}, skipping!".format(
                            os.path.split(png_file)[1], os.path.split(src_file)[1]))
                        found_preexisting = True
                        skipped_count += 1
                        pbar.update()
                        continue
                plot_xcorr_file_clock_analysis(src_file, dataset, time_window, snr_threshold, pearson_cutoff_factor,
                                               png_file=png_file, show=False, underlay_rcf_xcorr=underlay_rcf_xcorr,
                                               title_tag=title_tag, settings=settings)
            else:
                plot_xcorr_file_clock_analysis(src_file, dataset, time_window, snr_threshold, pearson_cutoff_factor,
                                               underlay_rcf_xcorr=underlay_rcf_xcorr, title_tag=title_tag,
                                               settings=settings)
            success_count += 1
            pbar.update()

        except Exception as e:
            tqdm.write("ERROR processing file {}".format(src_file))
            failed_files.append((src_file, str(e)))

        # Python 2 does not handle circular references, so it helps to explicitly clean up.
        if PY2:
            gc.collect()

    pbar.close()

    if found_preexisting:
        print("Some files were skipped because pre-existing matching png files were up-to-date.\n"
              "Remove png files to force regeneration.")

    return failed_files


def batch_process_folder(folder_name, dataset, time_window, snr_threshold, pearson_cutoff_factor=0.5, save_plots=True):
    """
    Process all the .nc files in a given folder into graphical visualizations.

    :param folder_name: Path to process containing .nc files
    :type folder_name: str
    :param dataset: Dataset to be used to ascertain the distance between stations.
    :type dataset: FederatedASDFDataset
    :param time_window: Lag time window to plot (plus or minus this value in seconds)
    :type time_window: float
    :param snr_threshold: Minimum signal to noise ratio for samples to be included into the clock lag estimate
    :type snr_threshold: float
    :param save_plots: Whether to save plots to file, defaults to True
    :param save_plots: bool, optional
    """
    src_files = glob.glob(os.path.join(folder_name, '*.nc'))
    print("Found {} .nc files in {}".format(len(src_files), folder_name))

    failed_files = batch_process_xcorr(src_files, dataset, time_window=time_window, snr_threshold=snr_threshold,
                                       pearson_cutoff_factor=pearson_cutoff_factor, save_plots=save_plots)
    _report_failed_files(failed_files)


def _report_failed_files(failed_files):
    if failed_files:
        print("The following files experienced errors:")
        for fname, err_msg in failed_files:
            print(" File: " + fname)
            if err_msg:
                print("Error: " + err_msg)


@click.command()
@click.argument('paths', type=click.Path('r'), required=True, nargs=-1)
@click.option('--dataset', default=None, type=click.Path('r'), help='Federated ASDF dataset filename')
@click.option('--time-window', default=300, type=int, show_default=True, help='Duration of time lag window to consider')
@click.option('--snr-threshold', default=6, type=float, show_default=True,
              help='Minimum sample SNR to include in clock correction estimate')
def main(paths, dataset, time_window, snr_threshold):
    """
    Main entry point for running clock analysis and CCF visualization on batch of station
    cross-correlation results.

    :param paths: Space separated list of files and folders to run analysis on. Files should be .nc files.
    :type paths: str or path
    :param dataset: Path to Federated ASDF dataset specification file
    :type dataset: str or path
    :param time_window: Time lag window size (x-axis extent) in seconds to use for analysis
    :type time_window: int
    :param snr_threshold: Minimum sample SNR to include in clock correction estimate
    :type snr_threshold: float
    """

    # Hardwired default path for now for GA application.
    if dataset is None:
        dataset = "/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt"

    if not os.path.exists(dataset):
        print("ERROR: Dataset path {} not found!".format(dataset))
        exit(1)

    files = []
    dirs = []
    for p in paths:
        if os.path.isdir(p):
            dirs.append(p)
        elif os.path.isfile(p):
            files.append(p)
        else:
            # click argument checking should ensure this never happens
            print("ERROR: Path {} not found!".format(p))
            sys.exit(1)

    ds = FederatedASDFDataSet.FederatedASDFDataSet(dataset)
    failed_files = batch_process_xcorr(files, ds, time_window, snr_threshold)
    _report_failed_files(failed_files)
    for d in dirs:
        batch_process_folder(d, ds, time_window, snr_threshold)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
