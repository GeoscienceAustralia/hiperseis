#!/usr/bin/env python
"""Utility plotting functions for consistent and convenient plotting of RFs.
"""

import logging
import os
import time
from collections import defaultdict

import numpy as np
import scipy.signal
from PyPDF2 import PdfFileMerger
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import signal
import rf

# pylint: disable=invalid-name, logging-format-interpolation

logging.basicConfig()

def revert_baz(rf_stream):
    def apply(st):
        for tr in st:
            if('orig_back_azimuth' in tr.stats.keys()):
                tr.stats.back_azimuth = tr.stats.orig_back_azimuth
            # end it
        # end for

        st.sort(['back_azimuth'])
        return st
    # end func

    if(isinstance(rf_stream, list)):
        result = [apply(st.copy()) for st in rf_stream]
    else:
        result = apply(rf_stream.copy())
    # end if

    return result
# end func

def plot_rf_psd(rf_stream, ax, time_window=(-10.0, 30.0), min_slope_ratio=-1):
    if(time_window): rf_stream = rf_stream.copy().slice2(*time_window, reftime='onset')

    all_trace_lens = np.array([len(tr) for tr in rf_stream])
    most_common_len, _ = stats.mode(all_trace_lens, axis=None)
    psd_stream = rf.RFStream([tr for tr in rf_stream if len(tr) == most_common_len])

    # plot psd
    psd_array = []
    fbins = None
    for trace in psd_stream:
        fbins, psd = signal.welch(trace.data, fs=trace.stats.sampling_rate,
                                  detrend='linear')

        psd_array.append(psd)
        ax.loglog(fbins, psd, alpha=0.05, c='m')
    # end for

    if (len(psd_array)):
        psd_array = np.array(psd_array)
        psd_mean = np.nanmean(psd_array, axis=0)

        ax.loglog(fbins, psd_mean, alpha=1, c='m', lw=2, label='Mean PSD')
        ax.set_xlabel('Freq. [Hz]')
        ax.set_ylabel('Power Spectral Density [arb. units]')
        ax.text(x=0.9, y=0.85, s='{} Traces'.format(len(psd_array)), transform=ax.transAxes)
        ax.legend()
        ax.grid()
    # end if
# end func

def plot_rf_stack(rf_stream, time_window=(-10.0, 30.0), trace_height=0.2, stack_height=0.8, save_file=None, **kwargs):
    """Wrapper function of rf.RFStream.plot_rf() to help do RF plotting with consistent formatting and layout.

    :param rf_stream: RFStream to plot
    :type rf_stream: rf.RFStream
    :param time_window: Time window to plot, defaults to (-10.0, 30.0)
    :type time_window: tuple, optional
    :param trace_height: Height of a single trace (reduce to cram RFs closer together), defaults to 0.2
    :type trace_height: float, optional
    :param stack_height: Height of mean (stacked) RF at top of plot, defaults to 0.8
    :type stack_height: float, optional
    :param save_file: File to save resulting image into, defaults to None
    :type save_file: str to valid file path, optional
    :return: Figure handle to the stack plot
    :rtype: matplotlib.figure.Figure
    """
    # Ensure traces are stackable by ignoring those that don't conform to the predominant data shape

    rf_stream = revert_baz(rf_stream)

    if(time_window): rf_stream = rf_stream.copy().slice2(*time_window, reftime='onset')

    all_trace_lens = np.array([len(tr) for tr in rf_stream])
    most_common_len, _ = stats.mode(all_trace_lens, axis=None)
    stackable_stream = rf.RFStream([tr for tr in rf_stream if len(tr) == most_common_len])
    num_stackable = len(stackable_stream)
    if num_stackable < len(rf_stream):
        logging.warning('Removed {} traces from RF plot to make it stackable!'.format(num_stackable))
    # end if

    fig = stackable_stream.plot_rf(fillcolors=('#000000', '#a0a0a0'), trim=time_window,
                                   trace_height=trace_height, stack_height=stack_height,
                                   fname=save_file, show_vlines=True, **kwargs)
    return fig
# end func


def plot_station_rf_overlays(db_station, title=None, time_range=None):
    """Plot translucent overlaid RF traces for all traces in each channel, and overplot
    the mean signal of all the traces per channel.

    :param db_station: Dictionary with list of traces per channel for a given station.
    :type db_station: dict({str, list(RFTrace)})
    :param title: Plot title, defaults to None
    :type title: str, optional
    :param time_range: Min and max times for the horizontal axis
    :type time_range: Pair of float
    :return: Mean trace signal per channel
    :rtype: list(numpy.array)
    """
    num_channels = 0
    for ch, traces in db_station.items():
        if traces:
            num_channels += 1

    plt.figure(figsize=(16, 8*num_channels))
    colors = ["#8080a040", "#80a08040", "#a0808040"]
    min_x = 1e+20
    max_x = -1e20

    signal_means = []
    for i, (ch, traces) in enumerate(db_station.items()):
        if not traces:
            continue
        col = colors[i]
        plt.subplot(num_channels, 1, i + 1)
        sta = traces[0].stats.station
        for j, tr in enumerate(traces):
            lead_time = tr.stats.onset - tr.stats.starttime
            times = tr.times()
            plt.plot(times - lead_time, tr.data, '--', color=col, linewidth=2)
            mask = (~np.isnan(tr.data) & ~np.isinf(tr.data))
            if j == 0:
                data_mean = np.zeros_like(tr.data)
                data_mean[mask] = tr.data[mask]  # pylint: disable=unsupported-assignment-operation
                counts = mask.astype(np.float)
            else:
                data_mean[mask] += tr.data[mask]
                counts += mask.astype(np.float)
            # end if
        # end for
        data_mean = data_mean/counts
        data_mean[(counts == 0)] = np.nan
        signal_means.append(data_mean)
        plt.plot(tr.times() - lead_time, data_mean, color="#202020", linewidth=2)
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude (normalized)')
        plt.grid(linestyle=':', color="#80808020")
        title_text = '.'.join([sta, ch])
        if title is not None:
            title_text += ' ' + title
        plt.title(title_text, fontsize=14)
        x_lims = plt.xlim()
        min_x = min(min_x, x_lims[0])
        max_x = max(max_x, x_lims[1])
    # end for
    if time_range is not None:
        min_x = time_range[0]
        max_x = time_range[1]

    for i in range(num_channels):
        subfig = plt.subplot(num_channels, 1, i + 1)
        subfig.set_xlim((min_x, max_x))
    # end for

    return signal_means


def plot_hk_stack(k_grid, h_grid, hk_stack, title=None, save_file=None, num=None, clip_negative=True,
                  depth_colour_range=(20, 70), stack_ylabel = 'Moho depth'):
    """Plot H-k stack using data generated by function seismic.receiver_fn.rf_stacking.computed_weighted_stack().

    :param k_grid: Grid of k-values
    :type k_grid: Two-dimensional numpy.array
    :param h_grid: Grid of H-values
    :type h_grid: Two-dimensional numpy.array
    :param hk_stack: Grid of stacked RF sample values produced by function
        seismic.receiver_fn.rf_stacking.computed_weighted_stack()
    :type hk_stack: Two-dimensional numpy.array
    :param title: Title to add to the plot, defaults to None
    :type title: str, optional
    :param save_file: File name in which to save the plot, defaults to None
    :type save_file: str or Path, optional
    :param num: Number of RFs used to produce the stack, defaults to None
    :type num: int, optional
    :param clip_negative: Clip negative stack regions to zero, defaults to True
    :type clip_negative: bool, optional
    :return: Handle to the figure created for the plot
    :rtype: matplotlib.figure.Figure
    """
    cmap = plt.cm.get_cmap('rainbow')
    fig = plt.figure(figsize=(16, 12))
    if clip_negative:
        hk_stack = hk_stack.copy()
        hk_stack[hk_stack < 0] = 0
    # end if

    assert(depth_colour_range[0] < depth_colour_range[1])
    cptmax = None
    extend = 'neither'
    if(depth_colour_range[0]>np.min(h_grid) or depth_colour_range[1]<np.max(h_grid)):
        cptmax = np.nanmax(hk_stack[(h_grid>=depth_colour_range[0]) & \
                                    (h_grid<=depth_colour_range[1])])
        cmap.set_over('magenta')
        extend = 'max'
    else:
        cptmax = np.nanmax(hk_stack)
    # end if
    cs = plt.contourf(k_grid, h_grid, hk_stack, levels=20, cmap=cmap, vmin=0, vmax=cptmax,
                      extend=extend)
    for c in cs.collections:
        c.set_rasterized(True)
    # end for

    cb = plt.colorbar(cs)

    cb.ax.set_ylabel('Stack sum')
    #plt.contour(k_grid, h_grid, hk_stack, levels=10, colors='k', linewidths=1)
    plt.xlabel(r'$\kappa = \frac{V_p}{V_s}$ (ratio)', fontsize=14)
    plt.ylabel('H = %s (km)'%(stack_ylabel), fontsize=14)
    if title is not None:
        plt.title(title, fontsize=16)
    plt.tick_params(right=True, labelright=True, which='both')
    plt.tick_params(top=True, labeltop=True, which='both')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.minorticks_on()
    plt.xlim(np.min(k_grid), np.max(k_grid))
    plt.ylim(np.min(h_grid), np.max(h_grid))

    if num is not None:
        xl = plt.xlim()
        yl = plt.ylim()
        txt_x = xl[0] + 0.95*(xl[1] - xl[0])
        txt_y = yl[0] + 0.95*(yl[1] - yl[0])
        plt.text(txt_x, txt_y, "N = {}".format(num), horizontalalignment='right',
                 color="#ffffff", fontsize=16, fontweight='bold', rasterized=True)
    # end if

    if save_file is not None:
        tries = 10
        while tries > 0:
            try:
                tries -= 1
                plt.savefig(save_file, dpi=300)
                break
            except PermissionError:
                time.sleep(1)
                if tries == 0:
                    print("WARNING: Failed to save file {} due to permissions!".format(save_file))
                    break
            # end try
        # end while
    # end if

    return fig
# end func


def plot_rf_wheel(rf_stream, max_time=15.0, deg_per_unit_amplitude=45.0, plt_col='C0', title='',
                  figsize=(10, 10), cluster=False, cluster_col='#ff4000', layout=None, fontscaling=1.0):
    """Plot receiver functions around a polar plot with source direction used to position radial RF plot.

    :param rf_stream: Collection of RFs to plot. If passed as a list, then each stream in the list
        will be plotted on separate polar axes.
    :type rf_stream: rf.RFStream or list(rf.RFStream)
    :param max_time: maximum time relative to onset, defaults to 25.0
    :type max_time: float, optional
    :param deg_per_unit_amplitude: Azimuthal scaling factor for RF amplitude, defaults to 20
    :type deg_per_unit_amplitude: float, optional
    :param plt_col: Plot color for line and positive signal areas, defaults to 'C0'
    :type plt_col: str, optional
    :param title: Title for the overall plot, defaults to ''
    :type title: str, optional
    :param figsize: Size of figure area, defaults to (12, 12)
    :type figsize: tuple, optional
    :param cluster: Whether to add overlaid mean RF where there are many RFs close together.
    :type cluster: bool
    :param cluster_col: Color of clustered stacked overlay plots.
    :type cluster_col: matplotlib color specification
    :param layout: Arrangement of polar plots in grid. If None, then arranged in a column.
    :type layout: tuple(int, int)
    :return: Figure object
    :rtype: matplotlib.figure.Figure
    """

    rf_stream = revert_baz(rf_stream)

    if not isinstance(rf_stream, list):
        rf_stream = [rf_stream]
    if layout is None:
        layout = (len(rf_stream), 1)

    figsize = (figsize[0]*layout[1], figsize[1]*layout[0])
    fig = plt.figure(figsize=figsize)

    for n, stream in enumerate(rf_stream):
        if not stream:
            continue
        ax = plt.subplot(*tuple(list(layout) + [n + 1]), projection="polar")
        # Orient with north
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_rlabel_position(75)

        inner_radius = 0.4*max_time  # time units (e.g. sec)
        stream = stream.copy().trim2(0, max_time, reftime='onset')
        for i, tr in enumerate(stream):
            t = tr.times()
            rf_amp = tr.data
            back_azi = np.deg2rad(tr.stats.back_azimuth)
            azi_amp = back_azi - np.deg2rad(deg_per_unit_amplitude*rf_amp/
                                            np.linspace(1, (np.max(t) - np.min(t))/inner_radius, len(t)))
            plt.plot(azi_amp, t, color=plt_col, zorder=i+1)
            ax.fill_betweenx(t, azi_amp, back_azi, where=((azi_amp - back_azi) < 0), lw=0., facecolor=plt_col,
                             alpha=0.7, zorder=i+1)
            ax.fill_betweenx(t, azi_amp, back_azi, where=((azi_amp - back_azi) >= 0), lw=0., facecolor='#a0a0a080',
                             zorder=i+1)
        # end for

        ax.set_rorigin(-inner_radius)
        ax.set_rlim(0, max_time)
        ax.tick_params(labelsize=14*fontscaling)

        stream_meta = stream[0].stats
        target_id = '.'.join([stream_meta.network, stream_meta.station, stream_meta.location, stream_meta.channel])
        ax.text(0.5, 0.5, target_id, fontsize=14*fontscaling, horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes)

        if cluster:
            try:
                from sklearn.cluster import DBSCAN

                back_azis = np.array([tr.stats.back_azimuth for tr in stream])
                clustering = DBSCAN(eps=1.0).fit_predict(back_azis.reshape(-1, 1))
                cluster_data = defaultdict(list)
                for i, cl in enumerate(clustering):
                    if cl == -1:
                        continue
                    cluster_data[cl].append(stream[i])
                # end for

                zplus = i + 1
                for i, cl in enumerate(cluster_data.values()):
                    # Have to assume same time samples for each RFTrace.
                    t = cl[0].times()
                    mean_azi = np.deg2rad(np.mean([tr.stats.back_azimuth for tr in cl]))
                    mean_amp = np.mean([tr.data for tr in cl], axis=0)

                    azi_amp = mean_azi - np.deg2rad(deg_per_unit_amplitude*mean_amp/
                                                    np.linspace(1, (np.max(t) - np.min(t))/inner_radius, len(t)))
                    plt.plot(azi_amp, t, color=cluster_col, zorder=i+zplus)
                    ax.fill_betweenx(t, azi_amp, mean_azi, where=((azi_amp - mean_azi) < 0), lw=0.,
                                     facecolor=cluster_col, zorder=i+zplus)
                    ax.fill_betweenx(t, azi_amp, mean_azi, where=((azi_amp - mean_azi) >= 0), lw=0.,
                                     facecolor='#a0a0a080', zorder=i+zplus)
                # end for
            except Exception as e:
                logging.error("Clustering RFs failed with error: {}".format(str(e)))
            # end try
        # end if

        if title:
            fig.suptitle(title, fontsize=20*fontscaling)
        # end if

    # end for

    return fig
# end func

def plot_iir_filter_response(filter_band_hz, sampling_rate_hz, corners):
    """Plot one-way bandpass filter response in the frequency domain. If filter is used as zero-phase,
    the attenuation will be twice what is computed here.

    :param filter_band_hz: Pair of frequencies corresponding to low cutoff and high cutoff freqs
    :type filter_band_hz: tuple(float) of length 2 (i.e. pair)
    :param sampling_rate_hz: The sampling rate in Hz
    :type sampling_rate_hz: float
    :param corners: The order of the filter
    :type corners: int
    :return: Figure object
    :rtype: matplotlib.figure.Figure
    """
    nyq_freq = sampling_rate_hz/2.0
    f_low = filter_band_hz[0]/nyq_freq
    f_high = filter_band_hz[1]/nyq_freq
    # Assuming code in obspy.signal.filter.bandpass uses this same iirfilter design function.
    z, p, k = scipy.signal.iirfilter(corners, [f_low, f_high], btype='band', ftype='butter', output='zpk')
    num_freqs = int(np.ceil(2*sampling_rate_hz/filter_band_hz[0]))
    w, h = scipy.signal.freqz_zpk(z, p, k, fs=sampling_rate_hz, worN=num_freqs)

    fig = plt.figure(figsize=(16, 9))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.set_title('Bandpass filter frequency response: band ({}, {})/{} Hz, order {}'
                  .format(filter_band_hz[0], filter_band_hz[1], sampling_rate_hz, corners), fontsize=18)
    ax1.plot(w, 10*np.log10(np.abs(h)), alpha=0.8, color='C0', linewidth=2)
    # ax1.set_xlim(0, 4*filter_band_hz[1])
    ax1.set_ylim(-100.0, 5)
    ax1.set_ylabel('Amplitude (dB)', color='C0', fontsize=16)
    ax1.tick_params(labelsize=16)
    ax1.set_xlabel('Frequency (Hz)', fontsize=16)
    ylims = plt.ylim()
    rect = patches.Rectangle((filter_band_hz[0], ylims[0]), filter_band_hz[1] - filter_band_hz[0],
                             ylims[1] - ylims[0], color='#80ff8040', zorder=0)
    ax1.add_patch(rect)
    plt.axvline(filter_band_hz[0], color='#80ff8080', linestyle='--')
    plt.axvline(filter_band_hz[1], color='#80ff8080', linestyle='--')
    plt.grid(linestyle=':', color="#80808080")

    ax2 = ax1.twinx()
    angles = np.unwrap(np.angle(h))
    ax2.plot(w, angles, alpha=0.8, color='C1', linestyle='--', linewidth=2)
    # ax2.set_xlim(0, 4*filter_band_hz[1])
    ax2.set_ylabel('Angle (rad)', color='C1', fontsize=16)
    ax2.tick_params(labelsize=16)

    # plt.axis('tight')

    return fig
# end func


def plot_iir_impulse_response(filter_band_hz, sampling_rate_hz, corners, zero_phase=False, N=1000, blip_period=1.0):
    """Plot bandpass filter response to standard waveforms in the time domain - impulse (delta function,
    step function, square wave pulse. By default filter is applied one-way. Set `zero_phase=True` to
    plot two-way filter response.

    :param filter_band_hz: Pair of frequencies corresponding to low cutoff and high cutoff freqs
    :type filter_band_hz: tuple(float) of length 2 (i.e. pair)
    :param sampling_rate_hz: The sampling rate in Hz
    :type sampling_rate_hz: float
    :param zero_phase: If True, plot two-way signal response (zero phase), otherwise plot
        one-way signal response.
    :type zero_phase: bool
    :param N: Number of samples in the input test signals.
    :type N: int
    :param blip_period: Period of the 'blip' test signal (square wave pulse).
    :type blip_period: float
    :return: Figure object
    :rtype: matplotlib.figure.Figure
    """
    nyq_freq = sampling_rate_hz/2.0
    f_low = filter_band_hz[0]/nyq_freq
    f_high = filter_band_hz[1]/nyq_freq
    # Assuming code in obspy.signal.filter.bandpass uses this same iirfilter design function.
    z, p, k = scipy.signal.iirfilter(corners, [f_low, f_high], btype='band', ftype='butter', output='zpk')
    sos = scipy.signal.zpk2sos(z, p, k)

    times = (np.arange(N) - (N//2))/sampling_rate_hz
    impulse = scipy.signal.unit_impulse(N, idx='mid')
    step = np.zeros_like(times)
    step[times >= 0] = 1
    # step -= 0.5
    blip = np.zeros_like(times)
    blip[(times >= -blip_period/2) & (times < 0)] = 1
    blip[(times >= 0) & (times <= blip_period/2)] = -1

    if zero_phase:
        ir = scipy.signal.sosfiltfilt(sos, impulse)
        sr = scipy.signal.sosfiltfilt(sos, step)
        br = scipy.signal.sosfiltfilt(sos, blip)
    else:
        ir = scipy.signal.sosfilt(sos, impulse)
        sr = scipy.signal.sosfilt(sos, step)
        br = scipy.signal.sosfilt(sos, blip)
    # end if

    yrange = 1.2

    fig = plt.figure(figsize=(16, 9))
    plt.subplot(3, 1, 1)
    plt.plot(times, impulse)
    plt.plot(times, ir, alpha=0.8)
    plt.title("Impulse response")
    plt.ylabel("Amplitude")
    plt.ylim(-yrange, yrange)
    plt.grid(linestyle=':', color="#80808080")
    plt.legend(['Input', 'Response'])

    plt.subplot(3, 1, 2)
    plt.plot(times, step)
    plt.plot(times, sr, alpha=0.8)
    plt.title("Step response")
    plt.ylabel("Amplitude")
    plt.ylim(-yrange, yrange)
    plt.grid(linestyle=':', color="#80808080")

    plt.subplot(3, 1, 3)
    plt.plot(times, blip)
    plt.plot(times, br, alpha=0.8)
    plt.title("Square pulse response (period = {} sec)".format(blip_period))
    plt.ylabel("Amplitude")
    plt.xlabel("Time (s)")
    plt.ylim(-yrange, yrange)
    plt.grid(linestyle=':', color="#80808080")

    direction = '(one way)' if not zero_phase else '(two way)'
    plt.suptitle("Reponse characteristics for filter band ({}, {})/{} Hz, order {} {}"
                 .format(filter_band_hz[0], filter_band_hz[1],
                         sampling_rate_hz, corners, direction))
    return fig
# end func

def pdf_merge(file_list, output_filename):
    merger = PdfFileMerger(strict=False)

    for pdffile in file_list:
        fn, _ = os.path.splitext(os.path.basename(pdffile))
        bookmark = '.'.join(fn.split('.')[1:])
        merger.append(pdffile, bookmark)
    # end for
    merger.write(output_filename)
    merger.close()
# end func
