#!/bin/env python
"""
Description:
    Generates a waveform-analytics report on raw data in mseed or asdf format

References:

CreationDate:   18/01/2024
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     07/02/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os, sys

is_windows = sys.platform.startswith('win')
is_osx = sys.platform.startswith('darwin')
is_linux = sys.platform.startswith('linux')

if(is_windows): # suppress warnings on windows
    import warnings
    warnings.filterwarnings('ignore', \
        message='Valid PROJ data directory not found')
# end if

import numpy as np
from obspy import UTCDateTime
import click
from typing import Callable
import matplotlib.pyplot as plt
import matplotlib
from collections import defaultdict
from obspy.core.inventory.response import Response
import matplotlib
from matplotlib import mlab
from obspy.signal.invsim import cosine_taper
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from obspy import read_inventory
from obspy.core import Stream
from seismic.misc import split_list
from seismic.ASDFdatabase.utils import MseedIndex
from seismic.ASDFdatabase.utils import MAX_DATE, MIN_DATE
from seismic.inventory.response import ResponseFactory
from tqdm import tqdm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as md
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime
import psutil
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocess
from multiprocess import Manager, freeze_support
from seismic.ASDFdatabase.analytics.sun import day_night_coverage, DAY_NIGHT_SECONDS
import uuid
from functools import partial
from glob import glob

if(is_windows | is_osx): matplotlib.use('TKAgg')
else: matplotlib.use('Agg')

# Sunrise/sunset times at Alice Springs are used for daytime and nighttime coverage analyses of
# waveform data. Coverage statistics computed as such will be off by +/- 1 hr at most for stations
# within continental Australia. This is an acceptable tradeoff for enhanced usability where users
# are not required to manually input station coordinates.
ALICE_SPRINGS_LON = 133.88
ALICE_SPRINGS_LAT = -23.70

class ProgressTracker(object):
    def __init__(self, manager: Manager):
        self.max = manager.Value('max', 0)
        self.i = manager.Value('i', 0)
        self.lock = manager.Lock()
    # end func

    def initialize(self, max_value, initial_value=0):
        with self.lock:
            self.max.value = max_value
            self.i.value = initial_value
            # end with
    # end func

    def increment(self):
        with self.lock:
            if (self.i.value < self.max.value):
                self.i.value += 1
            # end if
        # end with
    # end func

    def now(self):
        with self.lock:
            return self.i.value, self.max.value
        # end with
    # end func
# end class

class StationAnalytics():
    def __init__(self,
                 get_time_range_func: Callable[[str, str, str, str], tuple],
                 get_waveforms_func: Callable[[str, str, str, str, UTCDateTime, UTCDateTime],
                 Stream],
                 response: Response,
                 output_folder,
                 prog_tracker: ProgressTracker,
                 start_time: UTCDateTime = None,
                 end_time: UTCDateTime = None,
                 nproc=1):

        self.get_time_range_func = get_time_range_func
        self.get_waveforms_func = get_waveforms_func
        self.response = response
        self.output_folder = output_folder
        self.progress_tracker = prog_tracker
        self.start_time = start_time
        self.end_time = end_time
        self.nproc = nproc
        self.ppsd_list = None
        self.AMP_DB_MIN = -200
        self.AMP_DB_MAX = -50
        self.PER_MIN = 1e-2
        self.channel_min_st = defaultdict(lambda:UTCDateTime('2200-01-01'))
        self.channel_max_et = defaultdict(lambda:UTCDateTime(-1))
        # Red, Orange, Yellow, Green
        self.cmap = LinearSegmentedColormap.from_list('royg',
                                                      [(1, 0, 0),
                                                       (1, 0.5, 0),
                                                       (1, 1, 0),
                                                       (0, 1, 0)], N=5)
        self.cmap_norm = matplotlib.colors.Normalize(vmin=0, vmax=100)

        # fetch noise-models
        self.master_nm_periods, self.master_lnm = get_nlnm()
        _, self.master_hnm = get_nhnm()
        # flip noise-model values and restrict them to data range
        self.master_nm_periods = self.master_nm_periods[::-1]
        self.master_lnm = self.master_lnm[::-1]
        self.master_hnm = self.master_hnm[::-1]

        # create temp-folder
        self.temp_folder = os.path.join(self.output_folder, str(uuid.uuid4()))
        os.makedirs(self.temp_folder, exist_ok=True)
    # end func

    def analyse_data(self,
                     network: str,
                     station: str,
                     location: str,
                     channel: str,
                     sampling_rate: int):
        self.resp_amplitudes = None
        self.NFFT = 2 ** 16
        self.num_freqs = self.NFFT // 2 + 1
        self.freqs = np.fft.fftfreq(self.NFFT, 1 / sampling_rate)[:self.num_freqs]
        self.freqs = self.freqs[1:][::-1]
        self.freqs[0] *= -1
        self.periods = 1. / self.freqs
        self.sparse_periods = None

        # setup period binnings
        self._setup_period_bins()

        # create array for response amplitudes
        resp_amplitudes, _ = self.response.get_evalresp_response(t_samp=1 / sampling_rate,
                                                                 nfft=self.NFFT, output="VEL")
        resp_amplitudes = resp_amplitudes[1:]
        resp_amplitudes = resp_amplitudes[::-1]
        self.resp_amplitudes = np.absolute(resp_amplitudes * np.conjugate(resp_amplitudes))

        # collate timespans to be allocated to each parallel process
        st, et = self.get_time_range_func(network,
                                          station,
                                          location,
                                          channel)
        if(st < self.channel_min_st[channel]): self.channel_min_st[channel] = st
        if(et > self.channel_max_et[channel]): self.channel_max_et[channel] = et

        if (self.start_time and self.start_time > st): st = self.start_time
        if (self.end_time and self.end_time < et): et = self.end_time

        day_st = UTCDateTime(year=st.year, month=st.month, day=st.day)
        day_et = UTCDateTime(year=et.year, month=et.month, day=et.day) + DAY_NIGHT_SECONDS

        self.st_list = [cst for cst in np.arange(day_st, day_et, DAY_NIGHT_SECONDS)]
        self.et_list = [cst + DAY_NIGHT_SECONDS for cst in np.arange(day_st, day_et, DAY_NIGHT_SECONDS)]

        proc_st_list = split_list(self.st_list, self.nproc)
        proc_et_list = split_list(self.et_list, self.nproc)

        # initialize progress-tracker
        self.progress_tracker.initialize(len(self.st_list))

        # launch parallel computations
        if (1):
            p = Pool(ncpus=self.nproc)
            p.map(partial(self._generate_psds, network, station, location, channel, sampling_rate),
                  proc_st_list, proc_et_list)
        else:
            self._generate_psds(network, station, location, channel, sampling_rate,
                                self.st_list, self.et_list)
        # end if
    # end func

    def _setup_period_bins(self):
        period_smoothing_width_octaves = 1.0
        period_step_octaves = 0.5
        period_step_factor = 2 ** period_step_octaves
        period_smoothing_width_factor = \
            2 ** period_smoothing_width_octaves

        period_limits = [self.periods[0], self.periods[-1]]
        per_left = (period_limits[0] /
                    (period_smoothing_width_factor ** 0.5))
        per_right = per_left * period_smoothing_width_factor
        per_center = np.sqrt(per_left * per_right)

        # build up lists
        per_octaves_center = [per_center]
        # do this for the whole period range and append the values to our lists
        while per_center < period_limits[1]:
            # move left edge of smoothing bin further
            per_left *= period_step_factor
            # determine right edge of smoothing bin
            per_right = per_left * period_smoothing_width_factor
            # determine center period of smoothing/binning
            per_center = np.sqrt(per_left * per_right)
            # append to lists
            per_octaves_center.append(per_center)
        # wend

        # initialize periods used for subsampling raw PSD periods
        self.sparse_periods = np.array(per_octaves_center)

        # restrict nm to current set of periods
        hclip = self.master_nm_periods <= self.periods[-1]
        self.nm_periods = self.master_nm_periods[hclip]
        self.lnm = self.master_lnm[hclip]
        self.hnm = self.master_hnm[hclip]

        # initialize nearest-neighbour indices for mapping:
        # i) raw PSD periods to sparse-periods
        # ii) noise-model periods to sparse-periods
        self.dense_to_sparse_indices = np.fabs(self.sparse_periods[:, None] -
                                               self.periods[None, :]).argmin(axis=-1)
        self.nm_to_sparse_indices = np.fabs(self.sparse_periods[:, None] -
                                            self.nm_periods[None, :]).argmin(axis=-1)

        self.sparse_nm_overlap_indices = np.where((self.sparse_periods >= self.nm_periods[0]) &
                                                  (self.sparse_periods <= self.nm_periods[-1]))[0]
    # end func

    def _generate_psds(self, network, station, location, channel, sampling_rate,
                       start_time_list, end_time_list):
        def fft_taper(data):
            """
            Cosine taper, 10 percent at each end (like done by [McNamara2004]_).

            .. warning::
                Inplace operation, so data should be float.
            """
            data *= cosine_taper(len(data), 0.2)
            return data
        # end func

        for start_time, end_time in zip(start_time_list, end_time_list):
            # print(start_time, end_time)
            stream = self.get_waveforms_func(network,
                                             station,
                                             location,
                                             channel,
                                             start_time,
                                             end_time)

            stream_len = len(stream)
            if (stream_len == 0):
                self.progress_tracker.increment()
                continue
            else:
                sr = stream[0].stats.sampling_rate
                if (sr != sampling_rate):
                    print('Warning: discrepant sampling rate found. Expected {}, but found {} '
                          'in trace {} ({} - {}). Moving along..'.format(sampling_rate, sr,
                                                                         stream[0].get_id(),
                                                                         stream[0].stats.starttime,
                                                                         stream[0].stats.endtime))
                    self.progress_tracker.increment()
                    continue
                # end if
            # end if

            samples_processed = 0
            day_coverage_fraction = night_coverage_fraction = 0
            spec = None
            if (stream_len > 0):
                for i in np.arange(stream_len):
                    samples_processed += len(stream[i].data)
                    _spec, _ = mlab.psd(stream[i].data.astype('float32'), NFFT=self.NFFT,
                                        Fs=sampling_rate,
                                        detrend=mlab.detrend_linear, window=fft_taper,
                                        noverlap=0, sides='onesided',
                                        scale_by_freq=True)

                    if (spec is None):
                        spec = _spec
                    else:
                        spec += _spec
                    # end if
                # end for

                day_stream = stream.slice(endtime=end_time - 1e-5, nearest_sample=False)
                day_coverage_fraction, night_coverage_fraction = day_night_coverage(day_stream,
                                                                                    ALICE_SPRINGS_LON,
                                                                                    ALICE_SPRINGS_LAT)
                # print(start_time, day_coverage_fraction, night_coverage_fraction)
            # end for
            spec /= stream_len

            spec = spec[1:]
            spec = spec[::-1]

            # make omega with the same conventions as spec
            w = 2.0 * np.pi * self.freqs
            spec = (w ** 2) * spec / self.resp_amplitudes

            dtiny = np.finfo(0.0).tiny
            spec[spec < dtiny] = dtiny

            # comvert to dB
            spec = 10. * np.log10(spec)

            # compute coverage
            coverage_fraction = np.min([samples_processed / sampling_rate / (end_time - start_time), 1.])

            # compare spec to low- and high-noise-model
            sparse_spec = spec[self.dense_to_sparse_indices]

            sparse_spec_deviation_indices = np.where(
                (self.lnm[self.nm_to_sparse_indices] > sparse_spec) | \
                (self.hnm[self.nm_to_sparse_indices] < sparse_spec))[0]

            # Only include overlapping periods between nm and sparse_spec while calculating
            # spectral deviation
            overlap_deviation_indices = np.intersect1d(sparse_spec_deviation_indices,
                                                       self.sparse_nm_overlap_indices)
            deviation_fraction = overlap_deviation_indices.shape[0] /\
                                 self.sparse_nm_overlap_indices.shape[0]

            ############################################
            # Plot results and write output npz
            ############################################
            nslc = '.'.join((network, station, location, channel))
            output_fn_stem = '{}.{}.{}'.format(nslc,
                                               start_time,
                                               end_time)
            output_fn_stem = output_fn_stem.replace(':', '__')
            output_fn_png = os.path.join(self.temp_folder, output_fn_stem + '.png')
            output_fn_npz = os.path.join(self.temp_folder, output_fn_stem + '.npz')

            # save to npz
            np.savez(output_fn_npz, bounds=[start_time.timestamp, end_time.timestamp],
                     coverage_fraction=coverage_fraction,
                     sparse_periods=self.sparse_periods,
                     sparse_spec=sparse_spec,
                     deviation_fraction=deviation_fraction,
                     day_coverage_fraction=day_coverage_fraction,
                     night_coverage_fraction=night_coverage_fraction)

            face_color = ((1 - deviation_fraction) + coverage_fraction) / 2.
            # plot results
            if (1):
                try:
                    fig, ax = plt.subplots()
                    fig.set_size_inches(2, 1)

                    ax.set_facecolor(self.cmap(face_color))

                    ax.semilogx(self.sparse_periods, sparse_spec, 'b', lw=1)
                    ax.semilogx(self.nm_periods, self.lnm, 'k', lw=1)
                    ax.semilogx(self.nm_periods, self.hnm, 'k', lw=1)

                    ax.set_ylim(self.AMP_DB_MIN, self.AMP_DB_MAX)
                    ax.set_xlim(self.PER_MIN, np.max(self.periods))
                    ax.grid(True, which="both", ls="-", lw=0.2)
                    ax.tick_params(labelsize=6)

                    ax.text(0.015, -180,
                            'Coverage: {:.2f} %'.format(
                                coverage_fraction * 100),
                            fontsize=4, ha='left')
                    ax.text(0.015, -192,
                            'Spectral conformity: {:.2f} %'.format(
                                (1 - deviation_fraction) * 100),
                            fontsize=4, ha='left')

                    ax.text(0.015, -85,
                            nslc,
                            fontsize=4, ha='left')
                    fig.suptitle(start_time.strftime("%Y-%m-%d"), fontsize=7,
                                 y=0.85)

                    fig.savefig(output_fn_png, dpi=150, bbox_inches="tight", pad_inches=0)
                    plt.close()
                except Exception as e:
                    print(e)
                # end try
            # end if

            self.progress_tracker.increment()
            cv, mv = self.progress_tracker.now()
            print('Processing data under {}: [{}/{} days] {:2.1f}%'.format(nslc, cv, mv,
                                                                           cv / mv * 100), end='\r')
            sys.stdout.flush()
        # end for

        return None
    # end func

    def process_results(self, channel, output_fn):
        mean_spec = None
        mean_deviation = None
        total_coverage = None
        spec_count = 0

        nper = namp = self.sparse_periods.shape[0]
        amp_bins = np.linspace(self.AMP_DB_MIN, self.AMP_DB_MAX, namp)
        ppsd = np.zeros((nper, namp))
        days = []
        day_coverage_fractions = []
        night_coverage_fractions = []
        # load png files and results
        per_ids = np.arange(nper)

        npz_files = glob(os.path.join(self.temp_folder, '*.npz'))
        npz_files_dict = defaultdict(list)
        png_dict = defaultdict(list)
        for npz_file in npz_files:
            toks = os.path.basename(npz_file).split('.')
            net, sta, loc, cha, start_time = toks[:5]

            if(cha != channel): continue # process only required channel
            start_time = UTCDateTime(start_time).timestamp
            npz_files_dict[start_time].append(npz_file)

            png_file = os.path.splitext(npz_file)[0] + '.png'
            if (os.path.exists(png_file)):
                png_dict[start_time] = plt.imread(png_file)
            # end if
        # end for

        for start_time in tqdm(sorted(npz_files_dict.keys()),
                               desc='Loading results: '):
            for npz_file in npz_files_dict[start_time]:
                results = np.load(npz_file)
                spec_count += 1

                amp_ids = np.digitize(results['sparse_spec'], amp_bins)
                amp_ids[amp_ids == namp] = namp - 1
                ppsd[per_ids, amp_ids] += 1

                if (mean_spec is None):
                    mean_spec = results['sparse_spec']
                    total_coverage = results['coverage_fraction']
                    mean_deviation = results['deviation_fraction']
                else:
                    mean_spec += results['sparse_spec']
                    total_coverage += results['coverage_fraction']
                    mean_deviation += results['deviation_fraction']
                # end if

                days.append(results['bounds'][0])
                day_coverage_fractions.append(results['day_coverage_fraction'])
                night_coverage_fractions.append(results['night_coverage_fraction'])
            # end for
        # end for

        days = np.array(days)
        day_coverage_fractions = np.array(day_coverage_fractions)
        night_coverage_fractions = np.array(night_coverage_fractions)

        # check if any data was processed at all
        if (total_coverage is None):
            print('Warning: No results found..')
            return
        # end if

        # compute mean spectrum and deviation over all days
        if (spec_count):
            mean_spec /= spec_count
            mean_deviation /= spec_count
        # end if

        # compute total coverage fraction
        total_coverage_fraction = total_coverage / float(spec_count)

        # compute health for summary plot
        health = (((1 - mean_deviation) + total_coverage_fraction) / 2.) * 100

        #########################
        # generate report
        #########################
        if (os.path.splitext(output_fn)[1].upper() != '.PDF'): output_fn += '.pdf'

        with PdfPages(output_fn) as pdf:
            if (1):
                # title and summary
                fig, axes = plt.subplots(3, 1)
                fig.set_size_inches(8, 11)

                ax1, ax2, ax3 = axes

                # remove axes
                ax1.get_xaxis().set_visible(False)
                ax1.get_yaxis().set_visible(False)
                # Set report heading and overall assessment
                title = 'Analytics Report for {}'.format(channel)

                ax1.text(0.5, 0.7, title, ha='center', va='center', fontsize=18)
                ax1.text(0.5, 0.6, 'Report generated on: {} (UTC)'. \
                         format(UTCDateTime.now().strftime("%Y-%m-%dT%H:%M:%S")),
                         ha='center', va='center', fontsize=8)

                ax1.set_facecolor(self.cmap(self.cmap_norm(health)))
                ax1.text(0.015, 0.45,
                         'Data Coverage: {:.2f} %'.format(
                             total_coverage_fraction * 100),
                         fontsize=12, ha='left', c='k')
                ax1.text(0.015, 0.35,
                         'Spectral Conformity: {:.2f} %'.format(
                             ((1 - mean_deviation) * 100)),
                         fontsize=12, ha='left', c='k')
                ax1.text(0.015, 0.25,
                         'Recording Interval: {} - {} ({:.2f} DAYS)'.format(\
                             UTCDateTime(self.channel_min_st[channel]).strftime("%Y-%m-%d"),
                             UTCDateTime(self.channel_max_et[channel]).strftime("%Y-%m-%d"),
                             (self.channel_max_et[channel] - self.channel_min_st[channel]) /
                             DAY_NIGHT_SECONDS),
                         fontsize=12, ha='left', c='k')
                ax1.text(0.015, 0.15,
                         """Data Processed: {:.2f} DAYS""".format(total_coverage),
                         fontsize=12, ha='left', c='k')

                ax1.text(0.015, 0.05,
                         """Health = {:.2f} %""".format(health),
                         fontsize=12, ha='left', c='k', weight='bold')

                ax1.text(0.3, 0.05,
                         '(data_cov. + spec_conformity)/2',
                         fontsize=12, ha='left', c='k')

                # health colorbar
                divider = make_axes_locatable(ax1)
                cbax = divider.append_axes('right', size='2.5%', pad=0.1)
                sm = plt.cm.ScalarMappable(cmap=self.cmap, norm=self.cmap_norm)
                sm.set_array([])
                cbar = fig.colorbar(sm, cax=cbax, orientation='vertical')
                cbar.ax.set_ylabel('Health [%]')

                # generate summary plot
                try:
                    # plot PPSD
                    ax2.semilogx(self.sparse_periods, mean_spec, 'b', lw=2, label='Mean PSD')

                    sm = ax2.pcolormesh(self.sparse_periods, amp_bins, ppsd.T, rasterized=True)

                    ax2.semilogx(self.nm_periods, self.lnm, 'k', lw=2,
                                 label='Noise Model')
                    ax2.semilogx(self.nm_periods, self.hnm, 'k', lw=2)

                    ax2.set_ylim(-200, -50)
                    ax2.set_xlim(1e-2, np.max(self.periods))
                    ax2.grid(True, which="both", ls="-", lw=0.2, c='grey')
                    ax2.tick_params(labelsize=10)
                    ax2.set_xlabel("Period [s]", fontsize=10)
                    ax2.set_ylabel("Amp. [m2/s4][dB]", fontsize=10)
                    ax2.legend(loc='upper right')
                    ax2.set_title('PPSD')

                    divider = make_axes_locatable(ax2)
                    cbax = divider.append_axes('right', size='2.5%', pad=0.1)

                    cbar = fig.colorbar(sm, cax=cbax, orientation='vertical')
                    cbar.ax.set_ylabel('[%]')

                    # plot daily data coverage over day/night
                    days = [datetime.datetime.fromtimestamp(ts) for ts in days]

                    ax3.stem(days, 100 * day_coverage_fractions, 'salmon', linefmt='salmon', markerfmt='None',
                             basefmt='none', label='Daytime')
                    ax3.stem(days, -100 * night_coverage_fractions, 'deepskyblue', linefmt='deepskyblue',
                             markerfmt='None', basefmt='none', label='Nighttime')

                    leg = ax3.legend(loc='lower center', bbox_to_anchor=(0.5, -0.4), ncol=2)
                    for ht, color in zip(leg.get_texts(), ['salmon', 'deepskyblue']): ht.set_color(color)
                    ax3.set_ylabel("Coverage [%]", fontsize=10)

                    xfmt = md.DateFormatter('%Y-%m-%d')
                    ax3.xaxis.set_major_formatter(xfmt)
                    ax3.tick_params(axis='x', labelrotation=45)
                    ax3.tick_params(axis="x", direction="in")
                    ax3.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: str(abs(x))))
                    ax3.set_title('Data Coverage')

                except Exception as e:
                    print(e)
                # end try

                fig.tight_layout()
                pdf.savefig(dpi=300, bbox_inches="tight")
                plt.close()
            # end if

            if (1):
                def add_image(ax, img):
                    im = OffsetImage(img, zoom=0.4)
                    im.image.axes = ax
                    ab = AnnotationBbox(im, (0, 0),
                                        xybox=(30, 0.0),
                                        frameon=False,
                                        xycoords='data',
                                        boxcoords="offset points",
                                        pad=0)
                    ax.add_artist(ab)
                # end func

                # generate grids of daily plots
                nrows = 11
                ncols = 4
                nplots = np.sum([len(img) > 0 for _, img in png_dict.items()])
                sorted_keys = sorted(list(png_dict.keys()))
                plots_per_page = nrows * ncols

                done = False
                from matplotlib.offsetbox import OffsetImage, AnnotationBbox

                for ipage in np.arange(0, nplots, plots_per_page):
                    fig, axes = plt.subplots(nrows, ncols)
                    fig.set_size_inches(8, 11)

                    for irow in np.arange(nrows):
                        for icol in np.arange(ncols):
                            axes[irow, icol].set_axis_off()
                        # end for
                    # end for

                    for irow in np.arange(nrows):
                        for icol in np.arange(ncols):
                            iplot = int(ipage + irow * ncols + icol)

                            if ((iplot) >= nplots):
                                done = True
                                break;
                            # end if

                            ax = axes[irow, icol]

                            img = png_dict[sorted_keys[iplot]]
                            if (len(img) > 0):
                                add_image(ax, img)
                            # end if
                        # end for
                        if (done): break
                    # end for
                    pdf.savefig(dpi=300, bbox_inches='tight')
                    plt.close()
                    if (done): break
                # end for
            # end if
        # end with
    # end func
# end class

def get_response(input_file, network=None, station=None, location=None, channel=None):
    result = None
    if('xml' in input_file.lower()):
        inv = None
        try:
            inv = read_inventory(input_file)
        except Exception as e:
            print('Failed to read inventory file {} with error: {}'.format(input_file, e))
        # end try

        if(inv is not None):
            inv = inv.select(network=network, station=station, location=location,
                             channel=channel)
            if(inv is not None):
                seedid = inv.get_contents()['channels'][0]
                resp_obj = inv.get_response(seedid,
                                            inv.networks[0].stations[0].channels[0].start_date)
                result = resp_obj
        # end if
    else:
        resp_name = 'resp'
        # read resp file
        resp_inv = None
        try:
            resp_inv = read_inventory(input_file, format='RESP')
        except Exception as e:
            print('Failed to read RESP file {} with error: {}'.format(input_file, e))
        # end try

        if(resp_inv is not None):
            rf = ResponseFactory()
            rf.CreateFromInventory(resp_name, resp_inv)

            result = rf.getResponse(resp_name)
        # end if
    # end if

    return result
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)
def groups():
  pass
# end func

def get_cpu_count(nproc):
    result = nproc
    cpu_count = psutil.cpu_count(logical=False)
    if (nproc > cpu_count or nproc == -1):
        result = cpu_count
    # end if
    return result
# end func

def select_channel(mi:MseedIndex, sd:UTCDateTime, ed:UTCDateTime)->list([]):
    meta_list = mi.get_stations(sd, ed)
    cov_dict = mi.get_channel_coverages()
    channel_meta_dict = defaultdict(list)
    channel_cov_dict = defaultdict(float)

    for meta in meta_list:
        cha = meta[3]
        nslc = '.'.join(meta[:4])
        channel_meta_dict[cha].append(meta)
        channel_cov_dict[cha] += cov_dict[nslc]
    # end for

    if(len(meta_list) == 0):
        raise RuntimeError('No stations found between {} -- {}. Aborting..'.format(sd, ed))
    elif(len(meta_list) == 1):
        return channel_meta_dict[list(channel_meta_dict.keys())[0]]
    else:
        answer_list = [str(i + 1) for i in np.arange(len(channel_meta_dict))]
        answer = None
        print('\n############################')
        print('# Multiple channels found: #')
        print('############################\n')
        sorted_channels = sorted(channel_meta_dict.keys())
        for i, cha in enumerate(sorted_channels):
            print('{}) {} ({:2.3f} days)'.format(answer_list[i],
                                                cha, channel_cov_dict[cha]/DAY_NIGHT_SECONDS))
        # end for

        while (answer not in answer_list):
            answer = input('Please select desired channel: ')
        # wend

        channel_meta_list = channel_meta_dict[sorted_channels[int(answer)-1]]

        print('\nData found under the following NSLC list will be analysed:\n')
        for channel_meta in channel_meta_list:
            nslc = '.'.join(channel_meta)
            print('{} ({:2.3f} days)'.format(nslc, cov_dict[nslc]/DAY_NIGHT_SECONDS))
        # end for
        print('\n')

        return channel_meta_list
    # end if
# end func

@click.command(name='mseed', context_settings=CONTEXT_SETTINGS)
@click.argument('mseed-folder', required=True,
                type=click.Path(exists=True))
@click.argument('mseed-pattern', required=True,
                type=str)
@click.argument('instrument-response', required=True,
                type=click.Path(exists=True))
@click.argument('output-folder', required=True,
                type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.option('--start-date', type=str, default=None, show_default=True,
              help="Start date in UTC format for processing data")
@click.option('--end-date', type=str, default=None, show_default=True,
              help="End date in UTC format for processing data")
@click.option('--nproc', type=int, default=-1, show_default=True,
              help="Number of parallel processes to use. Default is to use all available cores.")
def process_mseed(mseed_folder, mseed_pattern, instrument_response,
                  output_folder, start_date, end_date, nproc):
    """
    MSEED_FOLDER: Path to folder containing mseed files\n
    MSEED_PATTERN: File pattern to be used to capture files pertaining to specific channels.
                   Note that pattern must be specified within quotes. \n
    INSTRUMENT_RESPONSE: Path to instrument response in .resp format\n
    OUTPUT_FOLDER: Path to output folder\n
    """
    try:
        start_date = UTCDateTime(start_date) if start_date else None
        end_date   = UTCDateTime(end_date) if end_date else None
    except Exception as e:
        print(str(e))
        raise RuntimeError('Invalid start- or end-dates')
    # end try

    if(start_date and end_date and ((end_date.date - start_date.date).days<=0)):
        raise RuntimeError('Invalid start- and end-dates. Aborting..')
    # end if

    # instantiate MseedIndex
    print('Inspecting mseed files in {}..'.format(mseed_folder))
    mseed_index = MseedIndex(mseed_folder, mseed_pattern)

    # define bound getter methods
    def get_waveforms_func(net, sta, loc, cha, st, et):
        return mseed_index.get_waveforms(net, sta, loc, cha, st, et)
    # end func

    def get_time_range_func(net, sta, loc, cha):
        return mseed_index.get_time_range(net, sta, loc, cha)
    # end func

    sd = MIN_DATE if start_date is None else start_date
    ed = MAX_DATE if end_date is None else end_date

    nproc = get_cpu_count(nproc)

    print('Loading response..')
    resp = get_response(instrument_response)
    if(resp is not None):
        print('Found response: {}'.format(resp))
        sys.stdout.flush()
    else:
        raise(RuntimeError('No instrument response found. Aborting..'))
    # end if

    # Recordings often feature garbled network, station and location codes.
    # We only consider unique channel codes as meaningful and as such, get a list of
    # net, sta, loc, cha combinations that share the same channel code. We also
    # account for the possibility of each combination of net, sta, loc, cha to have
    # recordings under different sampling rates
    meta_list = select_channel(mseed_index, sd, ed)
    sr_dict = mseed_index.get_channel_sampling_rates()

    # instantiate progress tracker
    manager = Manager()
    prog_tracker = ProgressTracker(manager)
    sa = StationAnalytics(get_time_range_func, get_waveforms_func,
                          resp, output_folder, prog_tracker, sd, ed, nproc)
    cha = None
    for meta in meta_list:
        net, sta, loc, cha = meta

        nslc = '.'.join(meta)
        sampling_rate = sr_dict[nslc]

        sa.analyse_data(net, sta, loc, cha, sampling_rate)
    # end for

    if(cha is not None): sa.process_results(cha, '')

    print('\nDone..')
# end func

@click.command(name='asdf', context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('network', required=True,
                type=str)
@click.argument('station', required=True,
                type=str)
@click.argument('location', required=True,
                type=str)
@click.argument('channel', required=True,
                type=str)
@click.argument('instrument-response', required=True,
                type=click.Path(exists=True))
@click.argument('sampling-rate', required=True,
                type=int)
@click.argument('output-folder', required=True,
                type=click.Path(exists=True))
@click.option('--start-date', type=str, default=None, show_default=True,
              help="Start date in UTC format for processing data")
@click.option('--end-date', type=str, default=None, show_default=True,
              help="End date in UTC format for processing data")
def process_asdf(asdf_source, network, station, location, channel, instrument_response,
                 sampling_rate, output_folder, start_date, end_date):
    """
    ASDF_SOURCE: Path to text file containing paths to ASDF files\n
    NETWORK: network code
    STATION: station code
    CHANNEL: channel code
    INSTRUMENT_RESPONSE: Path to inventory containing instrument response in
                         StationXML or .resp format\n
    SAMPLING_RATE: Sampling rate used to record the mssed files
    OUTPUT_FOLDER: Path to output folder\n
    """
    # import FederatedASDFDataSet locally to limit dependencies
    from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

    try:
        start_date = UTCDateTime(start_date) if start_date else None
        end_date   = UTCDateTime(end_date) if end_date else None
    except Exception as e:
        print(str(e))
        raise RuntimeError('Invalid start- or end-dates')
    # end try

    if(start_date and end_date and ((end_date.date - start_date.date).days<=0)):
        raise RuntimeError('Invalid start- and end-dates. Aborting..')
    # end if

    # instantiate FederatedASDFDataSet
    fds = FederatedASDFDataSet(asdf_source)

    sd = MIN_DATE if start_date is None else start_date
    ed = MAX_DATE if end_date is None else end_date
    meta_list = fds.get_stations(sd, ed, network=network, station=station, channel=channel)

    nslc = '{}.{}.{}.{}'.format(network, station, location, channel)
    if(len(meta_list) == 0):
        raise RuntimeError('No data found for {} between {} -- {}. Aborting..'.format(nslc, sd, ed))
    else:
        meta = meta_list[0]
    # end if
    net, sta, loc, cha = meta[:4]

    print('Loading response..')
    resp = get_response(instrument_response, net, sta, loc, cha)
    if(resp is not None): print('Found response: {}'.format(resp))
    else: raise(RuntimeError('No instrument response found. Aborting..'))

    # instantiate progress tracker
    manager = Manager()
    prog_tracker = ProgressTracker(manager)

    def get_waveforms_func(net, sta, loc, cha, st, et):
        return fds.get_waveforms(net, sta, loc, cha, st, et)
    # end func

    def get_time_range_func(net, sta, loc, cha):
        return fds.get_global_time_range(net, sta, loc, cha)
    # end func

    sa = StationAnalytics(get_time_range_func, get_waveforms_func,
                          prog_tracker, net, sta, loc, cha, sampling_rate, resp,
                          output_folder, sd, ed, nproc=1)

    report_fn = os.path.join(output_folder, '.'.join(meta[:4]) + '.pdf')
    sa.process_results(report_fn)
    print('Done..')
# end func

groups.add_command(process_mseed)
groups.add_command(process_asdf)

if __name__ == "__main__":
    # add support for process-based multiprocessing for a Windows .exe
    if(is_windows): freeze_support()
    if(is_osx): multiprocess.set_start_method('spawn')

    groups()
# end func
