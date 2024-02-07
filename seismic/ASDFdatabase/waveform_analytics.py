#!/bin/env python
"""
Description:
    Generate QAQC report on raw data in mseed or asdf format

References:

CreationDate:   18/01/2024
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     07/02/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
import numpy as np
from obspy import UTCDateTime
import click
from typing import Callable
import matplotlib.pyplot as plt
import matplotlib
from collections import defaultdict
from obspy.core.inventory.response import Response
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
from multiprocessing import Manager
from matplotlib.backends.backend_pdf import PdfPages
from joblib import Parallel, delayed
import psutil

os.sched_setaffinity(0, range(64))

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
                 prog_tracker: ProgressTracker,
                 network: str,
                 station: str,
                 location: str,
                 channel: str,
                 sampling_rate: int,
                 response: Response,
                 output_folder,
                 start_time: UTCDateTime = None,
                 end_time: UTCDateTime = None,
                 nproc=1):

        self.get_time_range_func = get_time_range_func
        self.get_waveforms_func = get_waveforms_func
        self.progress_tracker = prog_tracker
        self.network = network
        self.station = station
        self.location = location
        self.channel = channel
        self.sampling_rate = sampling_rate
        self.start_time = start_time
        self.end_time = end_time
        self.response = response
        self.output_folder = output_folder
        self.nproc = nproc
        self.ppsd_list = None
        self.resp_amplitudes = None
        self.NFFT = 2 ** 16
        self.num_freqs = self.NFFT // 2 + 1
        self.freqs = np.fft.fftfreq(self.NFFT, 1 / self.sampling_rate)[:self.num_freqs]
        self.freqs = self.freqs[1:][::-1]
        self.freqs[0] *= -1
        self.periods = 1. / self.freqs
        self.sparse_periods = None
        self.sparse_period_indices = None
        self.nm_periods, self.lnm = get_nlnm()
        _, self.hnm = get_nhnm()
        # Red, Orange, Yellow, Green
        self.cmap = LinearSegmentedColormap.from_list('royg',
                                                      [(1, 0, 0),
                                                       (1, 0.5, 0),
                                                       (1, 1, 0),
                                                       (0, 1, 0)], N=5)
        self.cmap_norm = matplotlib.colors.Normalize(vmin=0, vmax=100)
        self.overlap_denominator = None

        # flip noise-model values and restrict them to data range
        self.nm_periods = self.nm_periods[::-1]
        self.lnm = self.lnm[::-1]
        self.hnm = self.hnm[::-1]

        hclip = self.nm_periods <= self.periods[-1]
        self.nm_periods = self.nm_periods[hclip]
        self.lnm = self.lnm[hclip]
        self.hnm = self.hnm[hclip]
        self.st_list = None
        self.et_list = None

        # setup period binnings
        self._setup_period_bins()

        # create array for response amplitudes
        resp_amplitudes, _ = self.response.get_evalresp_response(t_samp=1 / self.sampling_rate,
                                                                 nfft=self.NFFT, output="VEL")
        resp_amplitudes = resp_amplitudes[1:]
        resp_amplitudes = resp_amplitudes[::-1]
        self.resp_amplitudes = np.absolute(resp_amplitudes * np.conjugate(resp_amplitudes))

        # collate timespans to be allocated to each parallel process
        day_seconds = 86400
        st, et = self.get_time_range_func(self.network,
                                          self.station,
                                          self.location,
                                          self.channel)
        if (self.start_time and self.start_time > st): st = self.start_time
        if (self.end_time and self.end_time < et): et = self.end_time

        day_st = UTCDateTime(year=st.year, month=st.month, day=st.day)
        day_et = UTCDateTime(year=et.year, month=et.month, day=et.day) + day_seconds

        self.st_list = [cst for cst in np.arange(day_st, day_et, day_seconds)]
        self.et_list = [cst + day_seconds for cst in np.arange(day_st, day_et, day_seconds)]

        proc_st_list = split_list(self.st_list, self.nproc)
        proc_et_list = split_list(self.et_list, self.nproc)

        # initialize progress-tracker
        self.progress_tracker.initialize(len(self.st_list))

        # launch parallel computations
        if (1):
            Parallel(n_jobs=self.nproc) \
                (delayed(self._generate_psds)(cst_list, cet_list) \
                 for cst_list, cet_list in zip(proc_st_list, proc_et_list))
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

        # initialize nearest-neighbour indices for mapping:
        # i) raw PSD periods to sparse-periods
        # ii) sparse-periods to noise-model periods
        self.dense_to_sparse_indices = np.fabs(self.sparse_periods[:, None] -
                                               self.periods[None, :]).argmin(axis=-1)
        self.sparse_to_nm_indices = np.fabs(self.nm_periods[None, :] -
                                            self.sparse_periods[:, None]).argmin(axis=-1)
        self.overlap_denominator = self.sparse_periods[
            np.where((self.sparse_periods >= self.nm_periods[0]) & \
                     (self.sparse_periods <= self.nm_periods[-1]))].shape[0]

    # end func

    def _generate_psds(self, start_time_list, end_time_list):
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
            stream = self.get_waveforms_func(self.network,
                                             self.station,
                                             self.location,
                                             self.channel,
                                             start_time,
                                             end_time)

            stream_len = len(stream)
            if (stream_len == 0):
                self.progress_tracker.increment()
                continue
            else:
                sr = stream[0].stats.sampling_rate
                if(sr != self.sampling_rate):
                    print('Warning: discrepant sampling rate found. Expected {}, but found {} '
                          'in trace {} ({} - {}). Moving along..'.format(self.sampling_rate, sr,
                                                                         stream[0].get_id(),
                                                                         stream[0].stats.starttime,
                                                                         stream[0].stats.endtime))
                    self.progress_tracker.increment()
                    continue
                # end if
            # end if

            samples_processed = 0
            spec = None
            if (stream_len > 0):
                for i in np.arange(stream_len):
                    samples_processed += len(stream[i].data)
                    _spec, _ = mlab.psd(stream[i].data.astype('float32'), NFFT=self.NFFT,
                                        Fs=self.sampling_rate,
                                        detrend=mlab.detrend_linear, window=fft_taper,
                                        noverlap=0., sides='onesided',
                                        scale_by_freq=True)
                    if (spec is None):
                        spec = _spec
                    else:
                        spec += _spec
                    # end if
                # end for
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
            coverage_fraction = np.min([samples_processed / self.sampling_rate / (end_time - start_time), 1.])

            # compare spec to low- and high-noise-model
            sparse_spec = spec[self.dense_to_sparse_indices]

            sparse_spec_deviation = np.where((self.lnm[self.sparse_to_nm_indices] > sparse_spec) | \
                                             (self.hnm[self.sparse_to_nm_indices] < sparse_spec))[0]

            ############################################
            # Plot results and write output npz
            ############################################
            deviation_fraction = len(sparse_spec_deviation) / self.overlap_denominator

            output_fn_stem = '{}.{}.{}.{}.{}.{}'.format(self.network,
                                                        self.station,
                                                        self.location,
                                                        self.channel,
                                                        start_time,
                                                        end_time)
            output_fn_png = os.path.join(self.output_folder, output_fn_stem + '.png')
            output_fn_npz = os.path.join(self.output_folder, output_fn_stem + '.npz')

            # save to npz
            np.savez_compressed(output_fn_npz, bounds=[start_time.timestamp, end_time.timestamp],
                                coverage_fraction=coverage_fraction,
                                sparse_periods=self.sparse_periods,
                                sparse_spec=sparse_spec,
                                deviation_fraction=deviation_fraction)

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

                    # ax.semilogx(self.sparse_periods[sparse_spec_deviation],
                    #            sparse_spec[sparse_spec_deviation], 'b',
                    #            linestyle=None, marker='+')

                    ax.set_ylim(-200, -50)
                    ax.set_xlim(1e-2, np.max(self.periods))
                    ax.grid(True, which="both", ls="-", lw=0.2)
                    ax.tick_params(labelsize=6)
                    # ax.set_xlabel("Period [s]", fontsize=5)
                    # ax.set_ylabel("Amp. [m2/s4][dB]", fontsize=5)
                    # ax.xaxis.set_label_coords(0.5, -0.04)
                    # ax.yaxis.set_label_coords(1.05, 0.5)

                    ax.text(0.015, -180,
                            'Coverage: {:.2f} %'.format(
                                coverage_fraction * 100),
                            fontsize=4, ha='left')
                    ax.text(0.015, -192,
                            'Spectral conformity: {:.2f} %'.format(
                                (1 - deviation_fraction) * 100),
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
            print('Progress: [{}/{}] {:2.1f}%'.format(cv, mv, cv / mv * 100), end='\r')
        # end for

        return None
    # end func

    def process_results(self, output_fn):
        png_dict = defaultdict(list)
        mean_spec = None
        mean_deviation = None
        total_coverage = None
        spec_count = 0

        # load png files and results
        for start_time, end_time in tqdm(zip(self.st_list, self.et_list),
                                         desc='Loading results: '):
            output_fn_stem = '{}.{}.{}.{}.{}.{}'.format(self.network,
                                                        self.station,
                                                        self.location,
                                                        self.channel,
                                                        start_time,
                                                        end_time)
            output_fn_png = os.path.join(self.output_folder, output_fn_stem + '.png')
            output_fn_npz = os.path.join(self.output_folder, output_fn_stem + '.npz')

            if (os.path.exists(output_fn_png)):
                png_dict[start_time.timestamp] = plt.imread(output_fn_png)
            # end if

            if (os.path.exists(output_fn_npz)):
                results = np.load(output_fn_npz)
                spec_count += 1

                if (mean_spec is None):
                    mean_spec = results['sparse_spec']
                    total_coverage = results['coverage_fraction']
                    mean_deviation = results['deviation_fraction']
                else:
                    mean_spec += results['sparse_spec']
                    total_coverage += results['coverage_fraction']
                    mean_deviation += results['deviation_fraction']
                # end if
            # end if
        # end for

        # check if any data was processed at all
        if(total_coverage is None):
            print('Warning: No results found..')
            return
        # end if

        # compute mean spectrum and deviation over all days
        if (spec_count):
            mean_spec /= spec_count
            mean_deviation /= spec_count
        # end if

        # compute total coverage fraction
        total_coverage_fraction = total_coverage * 86400 / \
                                  (self.et_list[-1] - self.st_list[0])

        # compute health for summary plot
        health = (((1 - mean_deviation) + total_coverage_fraction) / 2.) * 100

        #########################
        # generate report
        #########################
        if (os.path.splitext(output_fn)[1].upper() != '.PDF'): output_fn += '.pdf'

        with PdfPages(output_fn) as pdf:
            if (1):
                # title and summary
                fig, axes = plt.subplots(2, 1)
                fig.set_size_inches(8, 11)

                ax1, ax2 = axes
                ax1.set_axis_off()
                title = 'Analytics Report for {}.{}.{}.{}'.format(self.network,
                                                                  self.station,
                                                                  self.location,
                                                                  self.channel)
                ax1.text(0.5, 0.7, title, ha='center', va='center', fontsize=15)

                # generate summary plot
                try:
                    ax2.set_facecolor(self.cmap(self.cmap_norm(health)))
                    ax2.semilogx(self.sparse_periods, mean_spec, 'b', lw=2, label='PSD')
                    ax2.semilogx(self.nm_periods, self.lnm, 'k', lw=2,
                                 label='Noise Model')
                    ax2.semilogx(self.nm_periods, self.hnm, 'k', lw=2)

                    ax2.set_ylim(-200, -50)
                    ax2.set_xlim(1e-2, np.max(self.periods))
                    ax2.grid(True, which="both", ls="-")
                    ax2.tick_params(labelsize=10)
                    ax2.set_xlabel("Period [s]", fontsize=10)
                    ax2.set_ylabel("Amp. [m2/s4][dB]", fontsize=10)
                    ax2.legend(loc='upper right')

                    summary_text = \
                        """Summary PSD\nTime range: {} - {} ({:.2f} DAYS)\nData Processed: {:.2f} HOURS' or {:.2f} DAYS' worth in total""". \
                            format(self.st_list[0].strftime("%Y-%m-%d"),
                                   self.et_list[-1].strftime("%Y-%m-%d"),
                                   (self.et_list[-1] - self.st_list[0]) / 86400,
                                   total_coverage * 24, total_coverage)

                    ax2.text(0.015, -180,
                             'Coverage: {:.2f} %'.format(
                                 total_coverage_fraction * 100),
                             fontsize=15, ha='left')
                    ax2.text(0.015, -192,
                             'Spectral conformity: {:.2f} %'.format(
                                 ((1 - mean_deviation) * 100)),
                             fontsize=15, ha='left')

                    ax2.set_title(summary_text)

                    cbax = fig.add_axes([0.125, 0.05, 0.775, 0.015])
                    sm = plt.cm.ScalarMappable(cmap=self.cmap, norm=self.cmap_norm)
                    sm.set_array([])
                    cbar = fig.colorbar(sm, cax=cbax, orientation='horizontal')
                    cbar.ax.set_xlabel('Health: (coverage + spectral_conformity)/2')
                except Exception as e:
                    print(e)
                # end try

                ax1.set_axis_off()
                ax1.text(0.05, 0.05, 'Report generated on: {}'. \
                         format(UTCDateTime.now().strftime("%Y-%m-%dT%H:%M:%S")),
                         fontsize=7, ha='left')
                pdf.savefig(dpi=300, bbox_inches="tight")
                plt.close()
                # end if

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

            # end if

            # generate grids of daily plots
            nrows = 11
            ncols = 4
            nplots = (self.et_list[-1] - self.st_list[0]) / 86400
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

                        key = self.st_list[iplot].timestamp
                        img = png_dict[key]
                        if (img is not None):
                            add_image(ax, img)
                        # end if
                    # end for
                    if (done): break
                # end for
                pdf.savefig(dpi=300, bbox_inches='tight')
                plt.close()
                if (done): break
                # end for
        # end with
    # end func
# end class

def get_response(resp_file):
    resp_name = 'resp'
    # read resp file
    resp_inv = None
    try:
        resp_inv = read_inventory(resp_file, format='RESP')
    except Exception as e:
        print('Failed to read RESP file {} with error: {}'.format(resp_file, e))
    # end try

    rf = ResponseFactory()
    rf.CreateFromInventory(resp_name, resp_inv)

    return rf.getResponse(resp_name)
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)
def groups():
  pass
# end func

def get_cpu_count(nproc):
    result = nproc
    cpu_count = psutil.cpu_count()
    if (nproc > cpu_count or nproc == -1):
        result = cpu_count
    # end if
    return result
# end func

def select_channel(meta_list):
    answer_list = [str(i + 1) for i in np.arange(len(meta_list))]
    answer = None
    print('\n############################')
    print('# Multiple channels found: #')
    print('############################\n')
    for i in np.arange(len(meta_list)): print('{}) {}'.format(answer_list[i],
                                                              '.'.join(meta_list[i][:4])))
    while (answer not in answer_list):
        answer = input('Please select desired channel: ')
    # wend
    return meta_list[int(answer) - 1][:4]
# end func

@click.command(name='mseed', context_settings=CONTEXT_SETTINGS)
@click.argument('mseed-folder', required=True,
                type=click.Path(exists=True))
@click.argument('mseed-pattern', required=True,
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
@click.option('--nproc', type=int, default=-1, show_default=True,
              help="Number of parallel processes use. Default is to use all available cores.")
def process_mseed(mseed_folder, mseed_pattern, instrument_response,
                  sampling_rate, output_folder, start_date, end_date,
                  nproc):
    """
    MSEDD_FOLDER: Path to folder containing mseed files\n
    MSEED_PATTERN: File pattern to be used to capture files pertaining to specific channels.
                   Note that pattern must be specified within quotes. \n
    INSTRUMENT_RESPONSE: Path to instrument response in .resp format\n
    SAMPLING_RATE: Sampling rate used to record the mssed files
    OUTPUT_FOLDER: Path to output folder\n
    """
    try:
        start_date = UTCDateTime(start_date) if start_date else None
        end_date   = UTCDateTime(end_date) if end_date else None
    except Exception as e:
        print(str(e))
        raise RuntimeError('Invalid start- or end-dates')
    # end try

    print('Loading response..')
    resp = get_response(instrument_response)
    if(resp is not None): print('Found response: {}'.format(resp))
    else: raise(RuntimeError('No instrument response found. Aborting..'))

    # instantiate MseedIndex
    print('Inspecting mseed files in {}..'.format(mseed_folder))
    mseed_index = MseedIndex(mseed_folder, mseed_pattern)

    sd = MIN_DATE if start_date is None else start_date
    ed = MAX_DATE if end_date is None else end_date
    meta_list = mseed_index.get_stations(sd, ed)
    meta = None
    if(len(meta_list) == 0):
        raise RuntimeError('No stations found between {} -- {}. Aborting..'.format(sd, ed))
    elif(len(meta_list) > 1):
        meta = select_channel(meta_list)
    else:
        meta = meta_list[0]
    # end if
    net, sta, loc, cha = meta

    # instantiate progress tracker
    manager = Manager()
    prog_tracker = ProgressTracker(manager)

    def get_waveforms_func(net, sta, loc, cha, st, et):
        return mseed_index.get_waveforms(net, sta, loc, cha, st, et)
    # end func

    def get_time_range_func(net, sta, loc, cha):
        return mseed_index.get_time_range(net, sta, loc, cha)
    # end func

    nproc = get_cpu_count(nproc)
    sa = StationAnalytics(get_time_range_func, get_waveforms_func,
                          prog_tracker, net, sta, loc, cha, sampling_rate, resp,
                          output_folder, sd, ed, nproc)

    report_fn = os.path.join(output_folder, '.'.join(meta) + '.pdf')
    sa.process_results(report_fn)
    print('Done..')
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
    INSTRUMENT_RESPONSE: Path to instrument response in .resp format\n
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

    print('Loading response..')
    resp = get_response(instrument_response)
    if(resp is not None): print('Found response: {}'.format(resp))
    else: raise(RuntimeError('No instrument response found. Aborting..'))

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
    groups()
# end func
